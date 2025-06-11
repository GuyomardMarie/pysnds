import setuptools


VERSION = 0.1

setuptools.setup(
    name="pysnds",
    version=VERSION,
    author="Marie Guyomard",
    author_email="marie.guyomard83@gmail.com",
    description="A library to identify and interpret targeted information in the SNDS.",
    packages=setuptools.find_packages(),
    include_package_data=True,
    package_data={
        "pysnds": ["BC_medical_codes.json"],
    })
