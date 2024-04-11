import setuptools

with open("requirements.txt", "r") as f:
    install_requires = f.readlines()

setuptools.setup(
    name="optimorph",
    version="0.1.0",
    author="ronan-dupont",
    description="This is the OptiMorph version in 1D working in a cluster using SWAN or XBeach.",
    license="",
    long_description=open("README.md", "r").read(),
    classifiers=[
    ],
    include_package_data=True,
    install_requires=install_requires,
    python_requires='<=3.12',
)
