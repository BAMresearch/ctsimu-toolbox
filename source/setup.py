from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='ctsimu',
      version='1.4.4',
      description='CTSimU Software Toolbox',
      url='https://www.ctsimu.forschung.fau.de',
      author='David Plotzki',
      author_email='david.plotzki@bam.de',
      long_description=long_description,
      long_description_content_type="text/markdown",
      packages=['ctsimu', 'ctsimu.processing_pipeline'],
      python_requires='>=3.2',
      install_requires=['numpy', 'scipy'],
      include_package_data=True,
      classifiers=[
        "Programming Language :: Python :: 3",
        "License :: Other/Proprietary License",
        "Operating System :: OS Independent",
        "Natural Language :: English",
        "Topic :: Multimedia :: Graphics :: Graphics Conversion"
      ],
      zip_safe=True)