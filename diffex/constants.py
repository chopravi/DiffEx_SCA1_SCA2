import pandas as pd
import os

def package_path(**kwargs):
    """
    Return the path to the local installation
    of diffex
    """
    import diffex
    package_path = diffex.__file__
    package_path = os.path.dirname(package_path)
    return package_path

def source_path(**kwargs):
    """
    Return the path to the local installation
    source of diffex. (Includes 'data', 'notebooks',
    'bin', 'docs', etc.)
    """
    import diffex
    dirname = os.path.dirname(diffex.__file__)
    source_path = os.path.dirname(dirname)
    return source_path
    
package_path = package_path()
source_path = source_path()
analysis_path = os.path.join(source_path, 'notebooks')


IUPHAR_Channels = pd.read_csv(os.path.join(analysis_path, 'IUPHAR_Channels.csv'), encoding = "ANSI", engine='python')
IUPHAR_Channels_names = IUPHAR_Channels['Gene name']