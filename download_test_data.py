
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
    NETRC = None
    warnings.warn('Could not load ~/.netrc file.  Attempting to pull from DR17.')

HOST = 'data.sdss.org'
if NETRC is not None and HOST not in NETRC.hosts:
    NETRC = None
    warnings.warn('Host data.sdss.org is not defined in your ~/.netrc file.  '
                  'Attempting to pull from DR17.')

def download_file(remote_root, usr, passwd, local_root, file, overwrite=False):
    """
    Thanks to 
    https://stackoverflow.com/questions/37573483/progress-bar-while-download-file-over-http-with-requests/37573701
    """
    #Beware of how this is joined!
    url = f'{remote_root}{file}'
    ofile = os.path.join(local_root, file)

    if os.path.isfile(ofile):
        if overwrite:
            warnings.warn(f'Overwriting existing file: {ofile}')
            os.remove(ofile)
        else:
            warnings.warn('{ofile} exists.  To overwrite, set overwrite=True.')
            return

    print(f'Downloading: {url}')
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
    if NETRC is None:
        usr, acc, passwd = None, None, None
    else:
        usr, acc, passwd = NETRC.authenticators(HOST)

    sas_dir = 'dr17' if NETRC is None else 'mangawork'
    sas_root = f'https://{HOST}/sas/{sas_dir}/manga/spectro'

    version = drp_test_version
    files = remote_data_files()
    plates = [f.split('-')[1] for f in files]

    local_root = remote_data_file()
    if not os.path.isdir(local_root):
        os.makedirs(local_root)

    # Get the spectral data
    for plate, f in zip(plates, files):
        url_root = f'{sas_root}/redux/{drp_test_version}/{plate}/stack/'
        try:
            download_file(url_root, usr, passwd, local_root, f, overwrite=overwrite)
        except Exception as e:
            print(str(e))
            continue

    # Get the DRPComplete file
    f = f'drpcomplete_{drp_test_version}.fits'
    url_root = f'{sas_root}/analysis/{drp_test_version}/{dap_test_version}/common/'
    try:
        download_file(url_root, usr, passwd, local_root, f, overwrite=overwrite)
    except Exception as e:
        print(str(e))

    # Get the DRPall file
    f = f'drpall-{drp_test_version}.fits'
    url_root = f'{sas_root}/redux/{drp_test_version}/'
    try:
        download_file(url_root, usr, passwd, local_root, f, overwrite=overwrite)
    except Exception as e:
        print(str(e))


if __name__ == '__main__':
    main()


