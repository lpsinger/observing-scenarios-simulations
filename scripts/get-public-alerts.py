#!/usr/bin/env python
"""Create table of O3 public alerts from GraceDB."""
import warnings
import urllib.request

from astropy.table import Table
from astropy.utils.data import download_file
from gracedb_sdk import Client
from ligo.skymap import distance
from ligo.skymap.io import read_sky_map
from ligo.skymap.moc import uniq2pixarea
from ligo.skymap.postprocess.crossmatch import crossmatch
from ligo.skymap.util import progress_map
from lxml.etree import parse as parse_xml
import numpy as np
import requests.exceptions

client = Client(force_noauth=True)


def get_params_for_group(voevent_xml, name):
    elems = voevent_xml.findall(f".//Group[@type='{name}']/Param") or {}
    return {e.attrib['name']: float(e.attrib['value']) for e in elems}


def get_skymap(url):
    # Try to download the multiorder sky map, since it will be faster.
    try:
        new_url = url.replace('.fits.gz', '.multiorder.fits')
        filename = download_file(new_url, cache=True)
    except urllib.request.HTTPError:
        filename = download_file(url, cache=True)
    return read_sky_map(filename, moc=True)


def get_skymap_stats(skymap):
    crossmatch_result = crossmatch(skymap, contours=(0.9,), cosmology=True)
    area, = crossmatch_result.contour_areas
    vol, = crossmatch_result.contour_vols
    dist = distance.marginal_ppf(
        0.5, uniq2pixarea(skymap['UNIQ']) * skymap['PROBDENSITY'],
        skymap['DISTMU'], skymap['DISTSIGMA'], skymap['DISTNORM'])
    return {'area(90)': area, 'vol(90)': vol, 'distance': dist}


def get_info(superevent):
    superevent_id = superevent['superevent_id']
    result = {'superevent_id': superevent_id}
    api = client.superevents[superevent_id]

    # Scan VOEvents in reverse order (newest to oldest).
    voevents = api.voevents.get()[::-1]
    for voevent in voevents:
        try:
            voevent_xml = parse_xml(api.files[voevent['filename']].get())
        except requests.exceptions.HTTPError as e:
            # Some VOEvents cannot be found because the files in GraceDB were
            # not exposed to the public. Skip them.
            if e.response.status_code == 404:
                warnings.warn(f'HTTP Error 404 for f{voevent["filename"]}')
                continue
            else:
                raise
    
        # Get source classification and source properties.
        for group in ['Classification', 'Properties']:
            for key, value in get_params_for_group(voevent_xml, group).items():
                result.setdefault(key, value)

        # Look for a BAYESTAR sky map.
        elem = voevent_xml.find(".//Param[@name='skymap_fits']")
        if elem is not None and 'bayestar' in elem.attrib['value'].lower():
            skymap = get_skymap(elem.attrib['value'])
            for key, value in get_skymap_stats(skymap).items():
                result.setdefault(key, value)

        if all(key in result for key in ['Terrestrial', 'area(90)']):
            return result
    else:
        raise RuntimeError(f'Missing some information for {superevent_id}')


if __name__ == '__main__':
    # CBC events only
    superevents = (s for s in client.superevents.search(query='O3 ~ADVNO')
                    if s['preferred_event_data']['group'] == 'CBC')

    table = Table(rows=progress_map(get_info, superevents, jobs=None))

    # Add most likely source classification
    classifications = ['BNS', 'NSBH', 'BBH', 'MassGap']
    idx = np.argmax(table[classifications].columns.values(), axis=0)
    table['classification'] = np.asarray(classifications)[idx]

    # Reassign MassGap to NSBH if HasNS >= 0.5 else BBH
    table['classification'][(table['classification'] == 'MassGap') & (table['HasNS'] >= 0.5)] = 'NSBH'
    table['classification'][(table['classification'] == 'MassGap') & (table['HasNS'] < 0.5)] = 'BBH'

    # Put columns in a nicer order
    table = table['superevent_id', 'classification', 'distance',
                  'area(90)', 'vol(90)', 'HasNS', 'HasRemnant',
                  'BNS', 'NSBH', 'BBH', 'MassGap', 'Terrestrial']

    # Put rows in a nicer order
    table.sort('superevent_id')

    table.write('public-alerts.dat', format='ascii.tab', overwrite=True)
