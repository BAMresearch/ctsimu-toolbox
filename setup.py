import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(name='ctsimu',
      version='1.5.3',
      description='CTSimU Software Toolbox',
      url='https://bamresearch.github.io/ctsimu-toolbox/',
      author='David Plotzki',
      author_email='davidplotzki@gmail.com',
      long_description=long_description,
      long_description_content_type="text/markdown",
      packages=['ctsimu', 'ctsimu.processing', 'ctsimu.image_analysis', 'ctsimu.ctsimu_evaluations'],
      python_requires='>=3.8',
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