
import os
import warnings
import tqdm
import requests
import netrc

from IPython import embed

from mangadap.tests.util import remote_data_file, remote_data_files
from mangadap.tests.util import drp_test_version, dap_test_version

try:
    NETRC = netrc.netrc()
except Exception as e:
    raise FileNotFoundError('Could not load ~/.netrc file.') from e

HOST = 'data.sdss.org'
if HOST not in NETRC.hosts:
    raise ValueError('Host data.sdss.org is not defined in your ~/.netrc file.')


def download_file(remote_root, usr, passwd, local_root, file, overwrite=False):
    """
    Thanks to 
    https://stackoverflow.com/questions/37573483/progress-bar-while-download-file-over-http-with-requests/37573701
    """
    #Beware of how this is joined!
    url = '{0}{1}'.format(remote_root, file)
    ofile = os.path.join(local_root, file)

    if os.path.isfile(ofile):
        if overwrite:
            warnings.warn('Overwriting existing file: {0}'.format(ofile))
            os.remove(ofile)
        else:
            warnings.warn('{0} exists.  To overwrite, set overwrite=True.'.format(ofile))
            return

    print('Downloading: {0}'.format(url))
    # Streaming, so we can iterate over the response.
    r = requests.get(url, stream=True, auth=(usr, passwd))
    if r.status_code == 404:
        raise ValueError('Requested URL returned 404')
    # Total size in bytes.
    total_size = int(r.headers.get('content-length', 0))
    block_size = 1024 #1 Kibibyte
    t = tqdm.tqdm(total=total_size, unit='iB', unit_scale=True)
    with open(ofile, 'wb') as f:
        for data in r.iter_content(block_size):
            t.update(len(data))
            f.write(data)
    t.close()
    if total_size != 0 and t.n != total_size:
        raise ValueError('Downloaded file may be corrupted.')


def main():

    overwrite = False
    usr, acc, passwd = NETRC.authenticators(HOST)

    version = drp_test_version
    files = remote_data_files()
    plates = [f.split('-')[1] for f in files]

    local_root = remote_data_file()
    if not os.path.isdir(local_root):
        os.makedirs(local_root)

    # Get the spectral data
    for plate, f in zip(plates, files):
        url_root = 'https://{0}/sas/mangawork/manga/spectro/redux/{1}/{2}/stack/'.format(
                        HOST, drp_test_version, plate)
        try:
            download_file(url_root, usr, passwd, local_root, f, overwrite=overwrite)
        except Exception as e:
            print(str(e))
            continue

    # Get the DRPComplete file
    f = 'drpcomplete_{0}.fits'.format(drp_test_version)
    url_root = 'https://{0}/sas/mangawork/manga/spectro/analysis/{1}/{2}/common/'.format(
                    HOST, drp_test_version, dap_test_version)
    try:
        download_file(url_root, usr, passwd, local_root, f, overwrite=overwrite)
    except Exception as e:
        print(str(e))

    # Get the v3_0_1 drpcomplete until a v3_1_1 version is ready!
    f = 'drpcomplete_v3_0_1.fits'
    url_root = 'https://{0}/sas/mangawork/manga/spectro/analysis/{1}/{2}/common/'.format(
                    HOST, 'v3_0_1', '3.0.1')
    try:
        download_file(url_root, usr, passwd, local_root, f, overwrite=overwrite)
    except Exception as e:
        print(str(e))

    # Get the DRPall file
    f = 'drpall-{0}.fits'.format(drp_test_version)
    url_root = 'https://{0}/sas/mangawork/manga/spectro/redux/{1}/'.format(HOST, drp_test_version)
    try:
        download_file(url_root, usr, passwd, local_root, f, overwrite=overwrite)
    except Exception as e:
        print(str(e))


if __name__ == '__main__':
    main()


