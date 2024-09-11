from setuptools import setup, find_packages
import codecs
import os
from pathlib import Path

module_dir = Path(__file__).parent
with codecs.open(os.path.join(module_dir, "README.md"), encoding="utf-8") as fh:
    long_description = "\n" + fh.read()
#long_description = (module_dir / "README.md").read_text() #gets the long description
VERSION = '0.0.1'
DESCRIPTION = 'molecule screening'
LONG_DESCRIPTION = 'A package that will screen the best ligand for a given receptor with its grid'

# Setting up
setup(
    name="sakshi",
    version=VERSION,
    author="vijay",
    author_email="<vi.gade@ufl.edu",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=long_description,
    packages=find_packages(),
    #install_requires=['opencv-python', 'pyautogui', 'pyaudio'],
    keywords=['python', 'molecule', 'ligand', 'screening'],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Developers",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ]
)
