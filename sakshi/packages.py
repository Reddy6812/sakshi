import subprocess
import sys

def install_packages(packages):
    """
    Installs a list of packages using pip.
    
    Args:
        packages (list): A list of package names to install.
    """
    try:
        # Construct the pip install command
        command = ['pip', 'install'] + packages
        # Run the pip install command
        subprocess.check_call(command)
        print("Successfully installed packages:", " ".join(packages))
    except subprocess.CalledProcessError as e:
        print(f"Failed to install packages. Error: {e}")

if __name__ == "__main__":
    # List of packages to install
    packages = [
        'PyYAML',
        #'openbabel',
        'openbabel-wheel==3.1.1.7',
        'vina',
        'tqdm',
        'torch',
        'rdkit-pypi',
        'scikit-learn',
        'pandas',
        'rdkit'
    ]
    
    install_packages(packages)