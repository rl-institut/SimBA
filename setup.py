from setuptools import setup, find_packages

setup(
    name="simba",
    version="1.0.0",
    description="Simulation toolbox for Bus Applications",
    url="https://github.com/rl-institut/eBus-Toolbox",
    author="Reiner Lemoine Institut",
    author_email='info@rl-institut.de',
    license="MIT",
    packages=find_packages(),
    package_data={},
    install_requires=["dill", "matplotlib", "numpy", "pandas", "spice_ev"],
)
