from setuptools import setup, find_packages


setup(
    name="mdcv",
    version="0.0.1",
    packages=find_packages(),
    install_requires=[
#                    "mdciao",
                     ],
    dependency_links=[
        "https://github.com/gph82/mdciao/tarball/master",
    ]
)

