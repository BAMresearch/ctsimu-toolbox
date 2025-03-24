import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(name='ctsimu',
      version='1.6.0',
      description='CTSimU Software Toolbox',
      url='https://github.com/BAMresearch/ctsimu-toolbox/',
      author='David Plotzki',
      author_email='davidplotzki@gmail.com',
      long_description=long_description,
      long_description_content_type="text/markdown",
      packages=['ctsimu', 'ctsimu.processing', 'ctsimu.scenario', 'ctsimu.image_analysis', 'ctsimu.evaluation'],
      python_requires='>=3.9',
      include_package_data=True,
      classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
        "Natural Language :: English",
        "Topic :: Multimedia :: Graphics :: Graphics Conversion"
      ],
      zip_safe=True)