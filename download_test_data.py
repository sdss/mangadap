
import os
import warnings
import tqdm
import requests
import netrc

from mangadap.config import defaults
from mangadap.tests.util import remote_data_files

try:
    NETRC = netrc.netrc()
except Exception as e:
    raise FileNotFoundError('Could not load ~/.netrc file.') from e

HOST = 'data.sdss.org'
if HOST not in NETRC.hosts:
    raise ValueError('Host data.sdss.org is not defined in your ~/.netrc file.')


def download_file(remote_root, usr, passwd, local_root, file, overwrite=True):
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
            raise FileExistsError('File already exists.  To overwrite, set overwrite=True.')

    print('Downloading: {0}'.format(url))
    # Streaming, so we can iterate over the response.
    r = requests.get(url, stream=True, auth=(usr, passwd))
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

    usr, acc, passwd = NETRC.authenticators(HOST)

    version = 'MPL-9'
    files = remote_data_files()
    plates = [f.split('-')[1] for f in files]

    local_root = os.path.join(defaults.dap_data_root(), 'remote')

    for plate, f in zip(plates, files):
        url_root = 'https://{0}/sas/mangawork/manga/spectro/redux/{1}/{2}/stack/'.format(
                        HOST, version, plate)
        download_file(url_root, usr, passwd, local_root, f)


if __name__ == '__main__':
    main()


