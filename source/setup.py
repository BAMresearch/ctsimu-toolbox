from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='ctsimu',
      version='1.5',
      description='CTSimU Software Toolbox',
      url='https://www.ctsimu.forschung.fau.de',
      author='David Plotzki',
      author_email='davidplotzki@gmail.com',
      long_description=long_description,
      long_description_content_type="text/markdown",
      packages=['ctsimu', 'ctsimu.processing_pipeline', 'ctsimu.image_analysis', 'ctsimu.ctsimu_evaluations'],
      python_requires='>=3.9',
      install_requires=['numpy', 'scipy'],
      include_package_data=True,
      license="Apache 2.0",
      classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
        "Natural Language :: English",
        "Topic :: Multimedia :: Graphics :: Graphics Conversion"
      ],
      zip_safe=True)