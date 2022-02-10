# -*- coding: UTF-8 -*-
"""
This module provides classes for the virtual processing of images.

* `Image` reads, stores, writes and handles image data.
* `ImageFile` gathers information about an image file: file name, data type,
  byte order. It is used to instruct the `Image.read()` and `Image.save()`
  routines.
* `ImageStack` represents a stack of images in the file system. It can be used
  in combination with a processing pipeline (see `ctsimu.processing`).
* `ImageROI` defines a pixel region of interest in an image.

Images
------
To import a single image, you can specify its file name in the constructor
and then use the `Image.read()` function to import it into the internal memory.
It will be stored in `Image.px` as a float64 NumPy array. When writing an
image using `Image.save()`, you have to specify the data type for the new file.

    from ctsimu.image import Image
    
    myImage = Image("example.tif")
    myImage.read()
    
    # Mirror horizontally:
    myImage.flipHorizontal()
    
    myImage.save("example_mirrored.raw", dataType="float32")


RAW File Handling
-----------------
To read raw image data, its dimensions, data type, byte order and header size
must be specified:

    from ctsimu.image import Image

    myImage = Image("example_mirrored.raw")
    myImage.read(width=501,
                 height=501,
                 dataType="float32",
                 byteOrder="little",
                 fileHeaderSize=0)

    # Export as big endian, uint16:
    myImage.save("example_converted.raw",
                 dataType="uint16",
                 byteOrder="big")

"""

import numpy
import os    # File and path handling
import sys   # To get native byte order ('little' or 'big' endian?)
import math
import copy
from numpy.random import default_rng

# Scipy:
# 'ndimage' class for image processing
# 'optimize' class for intensity fit
# 'signal' class for drift analysis using FFT Convolution
from scipy import ndimage, optimize, stats, signal, fft

from .helpers import *
from .primitives import *   # Vectors and Polygons
from .tiffy import tiff

# pixelHalfDiagonal: longest distance a pixel center can have from a line
# while still touching the line with a corner point:
pixelHalfDiagonal = 1.0/math.sqrt(2.0)

def isTIFF(filename: str) -> bool:
    """Check if file name signifies a TIFF image."""
    if filename is not None:
        if(filename.casefold().endswith('.tif') or filename.casefold().endswith('.tiff')):
            return True
    
    return False

def createImageStack(stack):
    """ Return an ImageStack object, if string is given. """
    if isinstance(stack, ImageStack):
        return stack
    elif isinstance(stack, str):
        return ImageStack(stack)
    elif stack is None:
        return None
    else:
        raise Exception("Not a valid image file stack definition: {}".format(stack))

class ImageFile:
    """Fundamental image file properties used for input and output."""

    def __init__(self, filename=None, dataType=None, byteOrder=None, flipByteOrder=False):
        self.filename  = None
        self.dataType  = None
        self.byteOrder = None   # 'little' or 'big' endian
        self.flipByteOrder = False

        self.setFilename(filename)
        self.setDataType(dataType)
        self.setByteOrder(byteOrder)
        self.setFlipByteOrder(flipByteOrder)

    def setFilename(self, filename):
        self.filename = filename

    def getFilename(self) -> str:
        return self.filename

    def getFileBasename(self) -> str:
        return os.path.basename(self.filename)

    def getDataType(self) -> str:
        return self.dataType

    def getByteOrder(self) -> str:
        return self.byteOrder

    def doFlipByteOrder(self) -> bool:
        return self.flipByteOrder

    def setDataType(self, dataType: str):
        """ Set data type, either from numpy.dtype object or string. """
        if isinstance(dataType, numpy.dtype):
            self.dataType = dataType
        elif dataType is None:
            self.dataType = None
        elif isinstance(dataType, str):  # from string
            dt = numpy.dtype(dataType)
            self.setDataType(dt)
        else:
            raise Exception("{} is generally not a valid data type.".format(dataType))

    def setByteOrder(self, byteOrder: str):
        """ Set endianness, do sanity check before. """
        if byteOrder=='little' or byteOrder=='big' or byteOrder==None:
            self.byteOrder = byteOrder
        else:
            raise Exception("{} is not a valid byte order. Must be 'little' or 'big'.".format(byteOrder))

    def setFlipByteOrder(self, flipByteOrder: bool):
        self.flipByteOrder = flipByteOrder

    def isInt(self) -> bool:
        """ True if data type is supported int data type. """
        return numpy.issubdtype(self.dataType, numpy.integer)

    def isFloat(self) -> bool:
        """ True if data type is supported float data type. """
        return numpy.issubdtype(self.dataType, numpy.floating)

class ImageROI:
    """ Defines a region of interest: upper left and lower right corner. """

    def __init__(self, x0, y0, x1, y1):
        self.x0 = 0
        self.y0 = 0
        self.x1 = 0
        self.y1 = 0
        self.set(x0, y0, x1, y1)

    def __str__(self):
        return "({x0}, {y0}) -- ({x1}, {y1})".format(x0=self.x0, y0=self.y0, x1=self.x1, y1=self.y1)

    def set(self, x0, y0, x1, y1):
        if x1 < x0:
            x0, x1 = x1, x0

        if y1 < y0:
            y0, y1 = y1, y0

        self.x0 = int(x0)
        self.y0 = int(y0)
        self.x1 = int(x1)
        self.y1 = int(y1)

    def width(self):
        return self.x1 - self.x0

    def height(self):
        return self.y1 - self.y0

    def area(self):
        return self.width()*self.height()

    def grow(self, amount):
        amount = int(amount)
        self.set(self.x0-amount, self.y0-amount, self.x1+amount, self.y1+amount)


class Image:
    """ Stores pixel data, provides image processing routines. """

    def __init__(self, inputFile=None, outputFile=None):
        self.inputFile  = None   # type ImageFile or string
        self.outputFile = None   # type ImageFile or string
        self.px         = 0  # 2D numpy array that contains the pixel values.
        self.height     = 0  # Image height in px.
        self.width      = 0  # Image width in px.
        self.index      = 0  # Slice number in a 3D volume.

        self.rotation     = None
        self.flipHorz     = False
        self.flipVert     = False

        self.n_accumulations = 0   # Counts number of accumulated pictures for averaging (mean)
        self.boundingBoxX0   = 0   # After cropping: bounding box offset relative to original image.
        self.boundingBoxY0   = 0
        self.resolution      = 1   # After binning: new resolution relative to original image.

        self.setInputFile(inputFile)
        self.setOutputFile(outputFile)

    def __add__(self, other):
        if self.dimensionsMatch(other):
            result = copy.deepcopy(self)
            result.px += other.px
            return result
        else:
            raise Exception("Cannot add images of different dimensions.")

    def __sub__(self, other):
        if self.dimensionsMatch(other):
            result = copy.deepcopy(self)
            result.px -= other.px
            return result
        else:
            raise Exception("Cannot subtract images of different dimensions.")

    def __mul__(self, other):
        if self.dimensionsMatch(other):
            result = copy.deepcopy(self)
            result.px *= other.px
            return result
        else:
            raise Exception("Cannot multiply images of different dimensions.")

    def __truediv__(self, other):
        if self.dimensionsMatch(other):
            result = copy.deepcopy(self)
            result.px[numpy.nonzero(other.px)] /= other.px[numpy.nonzero(other.px)]
            result.px = numpy.where(other.px==0, 0, result.px)
            return result
        else:
            raise Exception("Cannot divide images of different dimensions.")

    def __floordiv__(self, other):
        if self.dimensionsMatch(other):
            result = copy.deepcopy(self)
            result.px[numpy.nonzero(other.px)] //= other.px[numpy.nonzero(other.px)]
            result = numpy.where(other.px==0, 0, result.px)
            return result
        else:
            raise Exception("Cannot divide images of different dimensions.")

    def __del__(self):
        """ Delete pixel map upon object destruction. """
        self.px =0

    def setInputFile(self, inputFile):
        """ Set input file properties from ImageFile object or string. """
        if isinstance(inputFile, ImageFile) or (inputFile is None):
            self.inputFile = inputFile
        elif isinstance(inputFile, str):  # string given
            self.inputFile = ImageFile(inputFile)
        else:
            raise Exception("{} is not a valid file identifier.")

    def setOutputFile(self, outputFile):
        """ Set output file properties from ImageFile object or string. """
        if isinstance(outputFile, ImageFile) or (outputFile is None):
            self.outputFile = outputFile
        elif isinstance(outputFile, str):  # string given
            self.outputFile = ImageFile(outputFile)
        else:
            raise Exception("{} is not a valid file identifier.")

    def setHeight(self, height):
        """ Set image height in px. """
        self.height = height

    def setWidth(self, width):
        """ Set image width in px. """
        self.width = width

    def setIndex(self, index):
        """ Set image index position in 3D stack (in px). """
        self.index = index

    def shape(self, width, height, index=0, dataType=None, value=0):
        """ Re-format image to given dimensions and data type. """
        self.setWidth(width)
        self.setHeight(height)
        self.setIndex(index)

        if dataType is None:
            dataType = self.getInternalDataType()

        self.erase(value=0, dataType=dataType)

    def shapeLike(self, otherImg, dataType=None):
        self.setWidth(otherImg.getWidth())
        self.setHeight(otherImg.getHeight())
        self.setIndex(otherImg.getIndex())

        if dataType is None:
            dataType = otherImg.getInternalDataType()

        self.erase(value=0, dataType=dataType)

    def erase(self, value=0, dataType=None):
        """ Set all pixels to 'value'. """
        w = self.getWidth()
        h = self.getHeight()

        if dataType is None:
            dataType = self.getInternalDataType()

        self.px = 0
        self.px = numpy.full((h, w), fill_value=value, dtype=dataType)
   
    def getPixelMap(self):
        return self.px

    def setPixelMap(self, px):
        self.px = px

    def setPixel(self, x, y, value):
        self.px[y][x] = value

    def getPixel(self, x, y):
        return self.px[y][x]

    def isSet(self):
        """ Check if image has a valid width and height. """
        if(self.getHeight() > 0):
            if(self.getWidth() > 0):
                return True

        return False

    def contains(self, x, y):
        """ Check if (x, y) is within image dimensions. """
        if x >= 0:
            if y >= 0:
                if x < self.getWidth():
                    if y < self.getHeight():
                        return True

        return False

    def getWidth(self):
        return self.width

    def getHeight(self):
        return self.height

    def getNPixels(self):
        """ Calculate number of pixels in image. """
        return (self.getWidth() * self.getHeight())

    def getIndex(self):
        return self.index

    def getBoundingBoxX0(self):
        return self.boundingBoxX0

    def getBoundingBoxY0(self):
        return self.boundingBoxY0

    def getResolution(self):
        return self.resolution

    def getFileByteOrder(self):
        return self.fileByteOrder

    def max(self, ROI=None):
        """ Return maximum intensity in image. """

        # Take full image if no ROI is given
        if ROI==None:
            return numpy.amax(self.px)

        return numpy.amax(self.px[ROI.y0:ROI.y1, ROI.x0:ROI.x1])

    def min(self, ROI=None):
        """ Return minimum intensity in image. """

        # Take full image if no ROI is given
        if ROI==None:
            return numpy.amin(self.px)

        return numpy.amin(self.px[ROI.y0:ROI.y1, ROI.x0:ROI.x1])

    def mean(self, ROI=None):
        """ Return arithmetic mean of the image grey values. """
        
        # Take full image if no ROI is given
        if ROI==None:
            return numpy.mean(self.px)

        return numpy.mean(self.px[ROI.y0:ROI.y1, ROI.x0:ROI.x1])

    def stdDev(self, ROI=None):
        """ Return the standard deviation of the image grey values. """

        # Take full image if no ROI is given
        if ROI==None:
            return numpy.std(self.px)

        return numpy.std(self.px[ROI.y0:ROI.y1, ROI.x0:ROI.x1])

    def centerOfMass(self):
        return ndimage.center_of_mass(self.px)

    def setRotation(self, rotation):
        self.rotation = rotation

    def getRotation(self):
        return self.rotation

    def rot90(self):
        if self.isSet():
            self.px = numpy.require(numpy.rot90(self.px, k=1), requirements=['C_CONTIGUOUS'])
            self.width, self.height = self.height, self.width

    def rot180(self):
        if self.isSet():
            self.px = numpy.require(numpy.rot90(self.px, k=2), requirements=['C_CONTIGUOUS'])

    def rot270(self):
        if self.isSet():
            self.px = numpy.require(numpy.rot90(self.px, k=-1), requirements=['C_CONTIGUOUS'])
            self.width, self.height = self.height, self.width

    def rotate(self, rotation):
        if rotation is None:
            rotation = self.rotation
        else:
            self.setRotation(rotation)

        if rotation == "90":
            self.rot90()
        elif rotation == "180":
            self.rot180()
        elif rotation == "270":
            self.rot270()

    def flipHorizontal(self):
        self.flipHorz = not self.flipHorz
        if self.isSet():
            self.px = numpy.require(numpy.fliplr(self.px), requirements=['C_CONTIGUOUS'])

    def flipVertical(self):
        self.flipVert = not self.flipVert
        if self.isSet():
            self.px = numpy.require(numpy.flipud(self.px), requirements=['C_CONTIGUOUS'])

    def setFlip(self, horz=False, vert=False):
        self.flipHorz = horz
        self.flipVert = vert

    def getHorizontalFlip(self):
        return self.flipHorz

    def getVerticalFlip(self):
        return self.flipVert

    def flip(self, horizontal=False, vertical=False):
        if horizontal:
            self.flipHorizontal()
        if vertical:
            self.flipVertical()

    def getInternalDataType(self):
        """ Data type used internally for all image data. """
        return numpy.dtype('float64')

    def containsPixelValue(self, value):
        """ Check if image contains a certain grey value. """
        return numpy.any(self.px == value)

    def dimensionsMatch(self, img):
        """ Check if image dimensions match with another image. """
        if self.isSet() and img.isSet():
            if(self.getHeight() == img.getHeight()):
                if(self.getWidth() == img.getWidth()):
                    return True

        raise Exception("Pixel dimensions do not match: {}x{} vs. {}x{}".format(self.getWidth(), self.getHeight(), img.getWidth(), img.getHeight()))
        
        return False

    def read(self, filename=None, width=None, height=None, index=0, dataType=None, byteOrder=None, fileHeaderSize=0, imageHeaderSize=0):
        """ Read TIFF or RAW, decide by file name. """
        if filename is None:
            filename = self.inputFile.getFilename()
        else:
            self.setInputFile(filename)

        # If no internal file name is specified, do nothing.
        if filename is None:
            return

        if isTIFF(self.inputFile.getFilename()):
            self.readTIFF(self.inputFile.doFlipByteOrder())
        else:
            self.readRAW(width=width, height=height, index=index, dataType=dataType, byteOrder=byteOrder, fileHeaderSize=fileHeaderSize, imageHeaderSize=imageHeaderSize)

    def readTIFF(self, flipByteOrder=False, obeyOrientation=True):
        """ Import TIFF file. """
        if os.path.isfile(self.inputFile.getFilename()):
            basename = self.inputFile.getFileBasename()
            
            tiffimg = tiff()
            tiffimg.read(self.inputFile.getFilename())
            img = tiffimg.imageData(subfile=0, channel=0, obeyOrientation=obeyOrientation)  # get a greyscale image from TIFF subfile 0
            width = tiffimg.getWidth(subfile=0)
            height = tiffimg.getHeight(subfile=0)

            self.inputFile.setDataType(img.dtype) 

            if flipByteOrder:
                img.byteswap(inplace=True)

            # Convert to internal data type for either int or float:
            self.px = img.astype(self.getInternalDataType())

            # Check if array in memory has the dimensions stated in the TIFF file:
            if((height == len(self.px)) and (width == len(self.px[0]))):
                self.setHeight(height)
                self.setWidth(width)
            else:
                raise Exception("Width ({}px) and height ({}px) from the TIFF header do not match the data width ({}px) and height ({}px) that has been read.".format(width, height, len(self.px[0]), len(self.px)))
        else:
            raise Exception("Can't find " + self.inputFile.getFilename())

    def readRAW(self, width, height, index=0, dataType=None, byteOrder=None, fileHeaderSize=0, imageHeaderSize=0):
        """ Import RAW image file. """
        if not isinstance(self.inputFile, ImageFile):
            raise Exception("No valid input file defined.")

        if dataType is None:
            dataType = self.inputFile.getDataType()
        else:
            self.inputFile.setDataType(dataType)

        if byteOrder is None:
            byteOrder = self.inputFile.getByteOrder()
            if byteOrder is None:
                byteOrder = sys.byteorder

        self.inputFile.setByteOrder(byteOrder)

        if os.path.isfile(self.inputFile.getFilename()):
            self.shape(width, height, index, self.inputFile.getDataType())

            basename = self.inputFile.getFileBasename()
            #log("Reading RAW file {}...".format(basename))

            byteOffset = fileHeaderSize + (index+1)*imageHeaderSize + index*(self.getNPixels() * self.inputFile.getDataType().itemsize)

            with open(self.inputFile.getFilename(), 'rb') as f:
                f.seek(byteOffset)
                self.px = numpy.fromfile(f, dtype=self.inputFile.getDataType(), count=self.getNPixels(), sep="")

            if len(self.px) > 0:
                # Treat endianness. If the native byte order of the system is different
                # than the given file byte order, the bytes are swapped in memory
                # so that it matches the native byte order.
                nativeEndian = sys.byteorder
                if nativeEndian == 'little':
                    if byteOrder == 'big':
                        self.px.byteswap(inplace=True)
                elif nativeEndian == 'big':
                    if byteOrder == 'little':
                        self.px.byteswap(inplace=True)

                # Convert to internal data type:
                self.px = self.px.astype(self.getInternalDataType())

                # Reshape to 2D array:
                self.px = numpy.reshape(self.px, (height, width))
            else:
                raise Exception("Error reading RAW file {f}.\nGot no data for index {idx}.".format(f=self.inputFile.getFilename(), idx=index))

        else:
            raise Exception("Can't find " + self.inputFile.getFilename())

    def getDataTypeClippingBoundaries(self, dataType):
        # Get clipping boundaries if grey values have to be
        # clipped to the interval supported by the int image type:
        clipMin = 0
        clipMax = 1
        if numpy.issubdtype(dataType, numpy.integer):
            intInfo   = numpy.iinfo(dataType)
            clipMin   = intInfo.min
            clipMax   = intInfo.max
        elif numpy.issubdtype(dataType, numpy.floating):
            floatInfo = numpy.finfo(dataType)
            clipMin   = floatInfo.min
            clipMax   = floatInfo.max

        return clipMin, clipMax

    def touchFolder(self, filename):
        """ Check if folder exists. Otherwise, create. """
        folder  = os.path.dirname(filename)
        if folder == "" or folder is None:
            folder = "."
        if not os.path.exists(folder):
            os.makedirs(folder)

    def save(self, filename=None, dataType=None, byteOrder=None, appendChunk=False, clipValues=True):
        """ Save image as TIFF or RAW. """
        if not isinstance(self.outputFile, ImageFile):
            self.outputFile = ImageFile()

        if (filename is None) or (filename == ""):
            filename = self.outputFile.getFilename()
            if (filename is None) or (filename == ""):
                raise Exception("No output file name specified.")
        else:
            self.outputFile.setFilename(filename)

        if dataType is None:
            dataType = self.outputFile.getDataType()
            if dataType is None:
                if isinstance(self.inputFile, ImageFile):
                    dataType = self.inputFile.getDataType()
                    if(dataType != None):
                        self.outputFile.setDataType(dataType)
                    else:
                        raise Exception("Please specify a data type for the output file: {filename}".format(filename=filename))
                else:
                    raise Exception("Please specify a data type for the output file: {filename}".format(filename=filename))
        else:
            self.outputFile.setDataType(dataType)

        if byteOrder is None:
            byteOrder = self.outputFile.getByteOrder()
            if byteOrder is None:
                if isinstance(self.inputFile, ImageFile):
                    byteOrder = self.inputFile.getByteOrder()
                    self.outputFile.setByteOrder(byteOrder)

            if byteOrder is None:
                byteOrder = "little"

        self.outputFile.setByteOrder(byteOrder)

        if isTIFF(filename):
            self.saveTIFF(filename, dataType, clipValues)
        else:
            self.saveRAW(filename, dataType, byteOrder, appendChunk, clipValues, addInfo=False)

    def saveTIFF(self, filename=None, dataType=None, clipValues=True):
        if (filename != None) and (len(filename) > 0):
            fileBaseName = os.path.basename(filename)
            if (fileBaseName == "") or (fileBaseName is None):
                raise Exception("No output file name specified for the image to be saved.")

            if dataType != None:
                if not isTIFF(filename):
                    filename += ".tif"

                self.touchFolder(filename)
                
                tiffdata = None
                if clipValues:  # Clipping
                    clipMin, clipMax = self.getDataTypeClippingBoundaries(dataType)
                    tiffdata = numpy.clip(self.px, clipMin, clipMax).astype(dataType)
                else:  # No clipping or float
                    tiffdata = self.px.astype(dataType)

                tiffimg = tiff()
                tiffimg.set(tiffdata)
                tiffimg.save(filename=filename, endian='little')
            else:
                raise Exception("Please specify a data type for the output file: {filename}".format(filename=filename))
        else:
            raise Exception("No output file name specified for the image to be saved.")
            
    def saveRAW(self, filename=None, dataType=None, byteOrder=None, appendChunk=False, clipValues=True, addInfo=False):
        if (filename != None) and (len(filename) > 0):
            fileBaseName = os.path.basename(filename)
            if (fileBaseName == "") or (fileBaseName is None):
                raise Exception("No output file name specified for the image to be saved.")

            if dataType != None:
                if byteOrder is None:
                    byteOrder = "little"

                # Reshape to 1D array and convert to file data type (from internal 64bit data type)
                outBytes = numpy.reshape(self.px, int(self.width)*int(self.height))

                if clipValues:  # Clipping
                    clipMin, clipMax = self.getDataTypeClippingBoundaries(dataType)
                    outBytes = numpy.clip(outBytes, clipMin, clipMax)

                outBytes = outBytes.astype(dataType)

                # Treat endianness. If the native byte order of the system is different
                # than the desired file byte order, the bytes are swapped in memory
                # before writing to disk.
                nativeEndian = sys.byteorder
                if nativeEndian == 'little':
                    if byteOrder  == 'big':
                        outBytes.byteswap(inplace=True)
                elif nativeEndian == 'big':
                    if byteOrder == 'little':
                        outBytes.byteswap(inplace=True)

                if addInfo:
                    shortEndian = "LE"
                    if byteOrder == "big":
                        shortEndian = "BE"

                    infoString = "_{width}x{height}_{dataType}_{endian}".format(width=self.width, height=self.height, dataType=dataType, endian=shortEndian)

                    basename, extension = os.path.splitext(filename)
                    filename = basename + infoString + extension

                self.touchFolder(filename)
                if not appendChunk:  # save as single raw file
                    with open(filename, 'w+b') as file:
                        file.write(outBytes)
                        file.close()
                    #outBytes.tofile(filename, sep="")
                else: # append to the bytes of the chunk file
                    with open(filename, 'a+b') as file:
                        file.write(outBytes)
                        file.close()
            else:
                raise Exception("Please specify a data type for the output file: {filename}".format(filename=filename))
        else:
            raise Exception("No output file name specified for the image to be saved.")

    def calcRelativeShift(self, referenceImage):
        if self.dimensionsMatch(referenceImage):
            # Convolution of this pixmap with the vertically and horizontally mirrored reference pixmap
            img1 = self.px - int(numpy.mean(self.px))
            img2 = referenceImage.getPixelMap() - numpy.mean(referenceImage.getPixelMap())

            convolution = signal.fftconvolve(img1, img2[::-1,::-1], mode='same')

            maximum = numpy.unravel_index(numpy.argmax(convolution), convolution.shape)

            return (maximum[1] - self.getWidth()/2, maximum[0] - self.getHeight()/2)
        else:
            raise Exception("Dimensions of image ({}, {}) and reference image ({}, {}) must match for convolution.".format(self.getWidth(), self.getHeight(), referenceImage.getWidth(), referenceImage.getHeight()))

    def getShiftedPixmap(self, xShift, yShift):
        return ndimage.interpolation.shift(self.px, (int(xShift), int(yShift)), mode='nearest')

    def accumulate(self, addImg, compensateShift=False, roiX0=None, roiY0=None, roiX1=None, roiY1=None):
        if (compensateShift == True) and (self.n_accumulations > 0):
            shift = (0, 0)

            if (roiX0 is None) or (roiY0 is None) or (roiX1 is None) or (roiY1 is None):
                shift = self.calcRelativeShift(addImg)
            else:
                # Crop image to drift ROI,
                croppedRef = copy.deepcopy(self)
                croppedRef.crop(x0=roiX0, y0=roiY0, x1=roiX1, y1=roiY1)

                croppedImg = copy.deepcopy(addImg)
                croppedImg.crop(x0=roiX0, y0=roiY0, x1=roiX1, y1=roiY1)

                shift = croppedImg.calcRelativeShift(croppedRef)

            log("Shift: {}".format(shift))
            shiftedPixMap = addImg.getShiftedPixmap(shift[1], shift[0])
            addImg.setPixelMap(shiftedPixMap)

        if self.n_accumulations == 0:
            self.setPixelMap(addImg.getPixelMap())
        else:
            if (self.dimensionsMatch(addImg)):
                self.px += addImg.getPixelMap()
            else:
                raise Exception("Current pixel dimensions ({currentX}x{currentY}) don't match dimensions of new file ({newX}x{newY}): {filename}".format(currentX=self.getWidth(), currentY=self.getHeight(), newX=addImg.getWidth(), newY=addImg.getHeight(), filename=addImg.inputFile.getFilename()))

        self.n_accumulations += 1

    def resetAccumulations(self):
        self.n_accumulations = 0

    def averageAccumulations(self):
        if self.n_accumulations > 1:
            self.px = self.px / self.n_accumulations
            log("Accumulated and averaged {} images.".format(self.n_accumulations))
            self.n_accumulations = 1

    def applyDark(self, dark):
        """ Apply dark image correction (offset). """
        if self.dimensionsMatch(dark):
            self.px = self.px - dark.getPixelMap()
        else:
            raise Exception("The dimensions of the image do not match the dimensions of the dark image for offset correction.")

    def applyFlatfield(self, ref, rescaleFactor=1):
        """ Apply flat field correction (free beam white image / gain correction). """
        if self.dimensionsMatch(ref):
            if(not ref.containsPixelValue(0)):  # avoid division by zero
                self.px = (self.px / ref.getPixelMap()) * float(rescaleFactor)
            else: # avoid division by zero
                self.px = (self.px / numpy.clip(ref.getPixelMap(), 0.1, None)) * float(rescaleFactor)
        else:
            raise Exception("The dimensions of the image do not match the dimensions of the flat image for flat field correction.")

    def verticalProfile(self, xPos):
        if xPos < self.getWidth():
            return numpy.ravel(self.px[:,xPos])
        else:
            raise Exception("Requested position for vertical profile is out of bounds: x={} in an image that has {} rows.".format(xPos, self.getWidth()))

    def verticalROIProfile(self, ROI):
        # Take full image if no ROI is given
        if ROI==None:
            ROI = ImageROI(0, 0, self.getWidth(), self.getHeight())

        slc = self.px[ROI.y0:ROI.y1, ROI.x0:ROI.x1]

        profile = slc.mean(axis=1)
        return numpy.ravel(profile)

    def horizontalProfile(self, yPos):
        if yPos < self.getHeight():
            return self.px[yPos]
        else:
            raise Exception("Requested position for horizontal profile is out of bounds: y={} in an image that has {} rows.".format(yPos, self.getHeight()))

    def horizontalROIProfile(self, ROI):
        # Take full image if no ROI is given
        if ROI==None:
            ROI = ImageROI(0, 0, self.getWidth(), self.getHeight())

        slc = self.px[ROI.y0:ROI.y1, ROI.x0:ROI.x1]

        profile = slc.mean(axis=0)
        return profile

    def pixelsInShape(self, shape, seedPoint=None, mode='center', calculateWeights=False):
        """ Returns all pixels in the given shape (of class Polygon). 

            mode:
              'center'   : a pixel's center must be within the shape to be accepted.
              'full'     : all corner points of a pixel must be within the shape to be accepted.
              'partial'  : only one corner point of a pixel must be within the shape to be accepted.

             calculateWeights:
               True      : includes weights in returned pixel coordinate tuples,
               False     : does not include weights in returned pixel coordinate tuples.
        """

        if seedPoint != None:
            seedX = int(round(seedPoint.x))
            seedY = int(round(seedPoint.y))
        else:
            # Start at point p1 of shape:
            seedX = int(shape.points[0].x)
            seedY = int(shape.points[0].y)

        # Make a map of visited pixels. A visited pixel will get value 1:
        visited = numpy.zeros_like(a=self.px, dtype=numpy.dtype('uint8'))

        # Collect all points that belong to the shape in a list:
        contributions = []

        stack = [] # stack of pixels to visit
        stack.append((seedX, seedY))

        # Add seed's neighors to the stack as well:
        for offsetX in [-1, 0, 1]:
            for offsetY in [-1, 0, 1]:
                if not (offsetX==0 and offsetY==0):
                    nx = seedX+offsetX
                    ny = seedY+offsetY
                    stack.append((nx, ny))

        while len(stack) > 0:
            pixel = stack.pop()
            x = pixel[0]
            y = pixel[1]

            if self.contains(x, y):
                if visited[y][x] == 0:
                    visited[y][x] = 1

                    # The pixel coordinate system is shifted by -0.5px against the shape coordinate system. Upper left pixel corner is its coordinate in the shape coordinate system.                   
                    inside = False

                    # Reserve names but set them up later only when they are needed.
                    center     = None
                    upperLeft  = None
                    upperRight = None
                    lowerLeft  = None
                    lowerRight = None

                    center = Vector(x+0.5, y+0.5, 0)

                    if mode == 'center':
                        inside = shape.isInside2D(center)
                    else:
                        upperLeft  = Vector(x,     y,     0)
                        upperRight = Vector(x+1,   y,     0)
                        lowerLeft  = Vector(x,     y+1,   0)
                        lowerRight = Vector(x+1,   y+1,   0)

                        if mode == 'full':
                            inside = shape.isInside2D(upperLeft) and shape.isInside2D(upperRight) and shape.isInside2D(lowerLeft) and shape.isInside2D(lowerRight)
                        elif mode == 'partial':
                            inside = True
                            calculateWeights = True
                    
                    if inside:
                        if calculateWeights:
                            # Calculate pixel weight from the area of the clipped pixel:
                            pixelPolygon = Polygon(upperLeft, upperRight, lowerRight, lowerLeft)  # Clockwise order because pixel CS is y-flipped.

                            clippedPixel = pixelPolygon.clip(shape)

                            weight = clippedPixel.area()

                            if weight > 0:
                                contributions.append((x, y, weight))
                            else:
                                continue
                        else:
                            contributions.append((x, y, 0))

                        # Now add neighbors to the stack:
                        for offsetX in [-1, 0, 1]:
                            for offsetY in [-1, 0, 1]:
                                if not (offsetX==0 and offsetY==0):
                                    nx = x+offsetX
                                    ny = y+offsetY
                                    stack.append((nx, ny))

        return contributions

    @staticmethod
    def getPixelWeight(x, y, clipPolygon):
        # Calculate pixel weight from the area of the clipped pixel:
        upperLeft  = Vector2D(x,   y)
        upperRight = Vector2D(x+1, y)
        lowerLeft  = Vector2D(x,   y+1)
        lowerRight = Vector2D(x+1, y+1)
        pixelPolygon = Polygon(upperLeft, upperRight, lowerRight, lowerLeft)  # Clockwise order because pixel CS is y-flipped.

        clippedPixel = pixelPolygon.clip(clipPolygon)
        weight = clippedPixel.area()

        return weight

    def meanGVinBin_polygonClipping(self, binCenter, sUnit, tUnit, sBoundary, tBoundary, binShape, weightFunction):
        """ Returns all pixels in the bin on the given vector s.

            binCenter:  center of bin in world CS
            s:   unit vector along profile axis
            t:   unit vector along width axis
        """

        roi_x0, roi_y0, roi_x1, roi_y1 = binShape.getBoundingBox()

        # Create a map with pixels' distances to the bin:
        # (measured parallel to s vector):
        roi_height = roi_y1 - roi_y0
        roi_width  = roi_x1 - roi_x0

        roi_xaxis = numpy.linspace(start=roi_x0, stop=roi_x1, num=roi_width+1, endpoint=True, dtype=numpy.dtype('float64'))
        roi_yaxis = numpy.linspace(start=roi_y0, stop=roi_y1, num=roi_height+1, endpoint=True, dtype=numpy.dtype('float64'))

        roi_gridx, roi_gridy = numpy.meshgrid(roi_xaxis, roi_yaxis)

        # Shift by half a pixel, because they must represent
        # pixel centers in shape coordinate system. Also,
        # origin should be the bin center:
        roi_gridx = roi_gridx + 0.5 - binCenter.x
        roi_gridy = roi_gridy + 0.5 - binCenter.y

        # Transform coordinates into bin coordinate system (s and t axes):
        bin_grid_dist_s = numpy.abs(roi_gridx*sUnit.x + roi_gridy*sUnit.y)
        #bin_grid_dist_t = numpy.abs(roi_gridx*tUnit.x + roi_gridy*tUnit.y)

        # Set those that are too far from bin center in s and t direction to zero:
        bin_grid_dist_s = numpy.where(bin_grid_dist_s < sBoundary, bin_grid_dist_s, 0)
        #bin_grid_dist_t = numpy.where(bin_grid_dist_t < tBoundary, bin_grid_dist_t, 0)
        #bin_grid_dist_mul = bin_grid_dist_s * bin_grid_dist_t
        #pixel_indices = numpy.nonzero(bin_grid_dist_mul)
        pixel_indices = numpy.nonzero(bin_grid_dist_s)
        pixels_x = pixel_indices[1] + roi_x0
        pixels_y = pixel_indices[0] + roi_y0

        weights = weightFunction(pixels_x, pixels_y, binShape)   # vectorized getPixelWeight()

        gvWeighted = self.px[pixels_y,pixels_x] * weights
        weightSum = numpy.sum(weights)
        meanGV = 0
        if weightSum > 0:
            meanGV = numpy.sum(gvWeighted) / weightSum

        return meanGV

    def meanGVinBin(self, binCenter, sUnit, tUnit, sBoundary, tBoundary, binShape, weightFunction):
        """ Returns all pixels in the bin on the given vector s.

            binCenter:  center of bin in world CS
            s:   unit vector along profile axis
            t:   unit vector along width axis
        """

        roi_x0, roi_y0, roi_x1, roi_y1 = binShape.getBoundingBox()

        # Create a map with pixels' distances to the bin:
        # (measured parallel to s vector):
        roi_height = roi_y1 - roi_y0
        roi_width  = roi_x1 - roi_x0

        roi_xaxis = numpy.linspace(start=roi_x0, stop=roi_x1, num=roi_width+1, endpoint=True, dtype=numpy.dtype('float64'))
        roi_yaxis = numpy.linspace(start=roi_y0, stop=roi_y1, num=roi_height+1, endpoint=True, dtype=numpy.dtype('float64'))

        roi_gridx, roi_gridy = numpy.meshgrid(roi_xaxis, roi_yaxis)

        # Shift by half a pixel, because they must represent
        # pixel centers in shape coordinate system. Also,
        # origin should be the bin center:
        roi_gridx = roi_gridx + 0.5 - binCenter.x
        roi_gridy = roi_gridy + 0.5 - binCenter.y

        # Transform coordinates into bin coordinate system (s and t axes):
        bin_grid_dist_s = numpy.abs(roi_gridx*sUnit.x + roi_gridy*sUnit.y)
        #bin_grid_dist_t = numpy.abs(roi_gridx*tUnit.x + roi_gridy*tUnit.y)

        # Set those that are too far from bin center in s and t direction to zero:
        #bin_grid_dist_s = numpy.where(bin_grid_dist_s < sBoundary, bin_grid_dist_s, 0)
        #bin_grid_dist_t = numpy.where(bin_grid_dist_t < tBoundary, bin_grid_dist_t, 0)
        #bin_grid_dist_mul = bin_grid_dist_s * bin_grid_dist_t
        #pixel_indices = numpy.nonzero(bin_grid_dist_mul)

        pixel_indices = numpy.nonzero(bin_grid_dist_s < sBoundary)
        weights = bin_grid_dist_s[pixel_indices]
        pixels_x = pixel_indices[1] + roi_x0
        pixels_y = pixel_indices[0] + roi_y0

        weights = weightFunction(pixels_x, pixels_y, binShape)   # vectorized getPixelWeight()

        gvWeighted = self.px[pixels_y,pixels_x] * weights
        weightSum = numpy.sum(weights)
        meanGV = 0
        if weightSum > 0:
            meanGV = numpy.sum(gvWeighted) / weightSum

        return meanGV

    """
    def lineProfile_projectPixelsIntoProfileBins(self, x0, y0, x1, y1, width=1, resolution=1):
        # Vector pointing in direction of the requested line:
        s = Vector(x1-x0+1, y1-y0+1, 0)   # +1 to fully include pixel (x1, y1)

        # Calculate vector t, perpendicular to s: t = s x z
        z = Vector(0, 0, 1)  
        t = s.cross(z)
        t.makeUnitVector()
        t.scale(0.5*width)

        # Define a rectangle along the line and its width, separated into two triangles.
        origin = Vector(x0, y0, 0)
        A = origin - t
        B = origin + s - t
        C = origin + s + t
        D = origin + t

        rect = Polygon(A, B, C, D)

        print("s: {}".format(s))
        print("t: {}".format(t))

        print(rect)

        ceilLength = math.ceil(s.length())

        nSamples = int( ceilLength / resolution ) + 1   # +1 for endpoint

        # Set a seed point at the center of the rectangle:
        t.scale(0.5)
        s.scale(0.5)
        seed = A + t + s

        # Make a list of unique pixel coordinates within this rectangle:
        pixelsInRect = self.pixelsInShape(shape=rect, seedPoint=seed)

        # Create a histogram:
        sPositions, sStepSize = numpy.linspace(start=0, stop=ceilLength, num=nSamples, endpoint=True, retstep=True)
        sCounts = numpy.zeros_like(a=sPositions, dtype=numpy.dtype('float64'))   # Number of contributions, for correct re-normalization, same datatype for efficiency during division later on...
        sSum = numpy.zeros_like(a=sPositions, dtype=numpy.dtype('float64'))      # The sum of all grey value contributions

        # Make s a unit vector to correctly calculate projections using the dot product:
        s.makeUnitVector()

       # print("shape of positions: {}".format(numpy.shape(sPositions)))

        print("{} pixels in rect.".format(len(pixelsInRect)))

        offset = Vector(0.5, 0.5, 0)

        for pixel in pixelsInRect:
            # Project this pixel onto the s vector (pointing in direction of the line):

            # Move to line origin:
            p = pixel - origin + offset

            # Position on s axis:
            sPos = p.dot(s)

            # Find bin where this grey value should be counted:
            binPos = int(math.floor(sPos / sStepSize))

            #print("({x}, {y}): sPos: {spos}, binPos: {binpos}".format(x=p.x, y=p.y, spos=sPos, binpos=binPos))

            sCounts[binPos] += 1
            sSum[binPos] += self.getPixel(int(pixel.x), int(pixel.y))

        # Replace zero counts by 1 to avoid div by zero:
        sCounts[sCounts==0] = 1

        sProfile = sSum / sCounts

        return sProfile, sPositions, sStepSize
    """

    def lineProfile(self, x0, y0, x1, y1, width=1, resolution=1):
        """ Find line profile by adding weighted contributions of pixel grey values
            into bins of size (width x resolution).

            We always work in the 'shape coordinate system' with its origin
            at (0, 0) in the upper left corner.
            Center of pixel (0, 0) has shape CS coordinates (0.5, 0.5).

            x0, y0, x1 and y1 are shape coordinates.

            Returned 'sPositions' array contains bin center positions.
            """

        # Vector pointing in direction of the requested line:
        s = Vector(x1-x0, y1-y0, 0)

        # Calculate vector t, perpendicular to s: t = s x z
        z = Vector(0, 0, 1)  
        t = s.cross(z)
        t.makeUnitVector()

        # Convert to 2D vectors:
        s = Vector2D(s.x, s.y)
        t = Vector2D(t.x, t.y)

        tUnit = copy.deepcopy(t)

        t.scale(0.5*width)  # t points from line origin half way in direction of width

        # Define a rectangle along the line and its width.
        origin = Vector2D(x0, y0)

        nSamples = math.ceil( s.length() / resolution ) #+ 1 # +1 for endpoint
        ceilLength = nSamples * resolution

        # Create a histogram:
        sPositions, sStepSize = numpy.linspace(start=0, stop=ceilLength, num=nSamples, endpoint=False, retstep=True)
        sProfile = numpy.zeros_like(a=sPositions, dtype=numpy.dtype('float64'))   # Grey value profile

        # Create a unit vector in s direction:
        sUnit = copy.deepcopy(s)
        sUnit.makeUnitVector()

        # Half a unit vector:
        binUnitHalf = copy.deepcopy(sUnit)
        binUnitHalf.scale(0.5*resolution)

        # Make s the length of a bin step (i.e. resolution unit)
        s.makeUnitVector()
        s.scale(resolution)

        rectPos = Vector2D(0, 0)

        # A pixel center can be this far from the binPos (bin center)
        # in s and t direction to still be accepted:
        sBoundary = (resolution/2) + pixelHalfDiagonal
        tBoundary = (width/2) + pixelHalfDiagonal

        # Vectorize the pixel weight function:
        weightFunction = numpy.vectorize(self.getPixelWeight, otypes=[numpy.float64])

        i = 0
        for b in range(nSamples):
            print("\rCalculating line profile... {:.1f}%".format(100.0*i/nSamples), end="")
            i += 1
            # Bin position on s axis:
            sPos = resolution*b

            # Construct a vector to the left point of the bin on the s axis:
            rectPos.setx(sUnit.x)
            rectPos.sety(sUnit.y)
            rectPos.scale(sPos)
            rectPos.add(origin)

            binPos = rectPos + binUnitHalf

            # Construct a rectangle that contains the area of this bin:
            A = rectPos - t
            B = rectPos + s - t
            C = rectPos + s + t
            D = rectPos + t

            binRect = Polygon(D, C, B, A)  # Clockwise order because pixel CS is y-flipped.

            # Get all pixels and their relative areas in this bin:
            #pixelsInBin = self.pixelsInShape(shape=binRect, seedPoint=rectPos, mode='partial', calculateWeights=True)

            meanGV = self.meanGVinBin(binCenter=binPos, sUnit=sUnit, tUnit=tUnit, sBoundary=sBoundary, tBoundary=tBoundary, binShape=binRect, weightFunction=weightFunction)

            sProfile[b] = meanGV

        # Shift the sPositions by half a bin size so that they represent bin centers:
        sPositions += 0.5*resolution

        print("\rCalculating line profile... 100%   ")
        return sProfile, sPositions, sStepSize
                
    def clip(self, lower, upper):
        """ Clip grey values to given boundary interval. """
        self.px = numpy.clip(self.px, lower, upper)

    def crop(self, x0, y0, x1, y1):
        """ Crop to given box (x0, y0)--(x1, y1). """
        if x0 > x1:
            x0,x1 = x1,x0

        if y0 > y1:
            y0,y1 = y1,y0

        if y1 > self.getHeight()  or  x1 > self.getWidth():
            raise Exception("Trying to crop beyond image boundaries.")

        self.boundingBoxX0 += x0
        self.boundingBoxY0 += y0

        self.px = self.px[int(y0):int(y1),int(x0):int(x1)]   # Array has shape [y][x]
        self.width  = int(x1 - x0)
        self.height = int(y1 - y0)

    def cropBorder(self, top=0, bottom=0, left=0, right=0):
        """ Crop away given border around image. """
        x0 = int(left)
        y0 = int(top)
        x1 = self.getWidth() - int(right)
        y1 = self.getHeight() - int(bottom)

        self.crop(x0, y0, x1, y1)

    def cropROIaroundPoint(self, centerX, centerY, roiWidth, roiHeight):
        """ Crop a region of interest, centerd around given point. """

        if roiWidth < 0:
            roiWidth = abs(roiWidth)
        if roiHeight < 0:
            roiHeight = abs(roiHeight)
        if roiWidth == 0 or roiHeight == 0:
            raise Exception("The region of interest should not be a square of size 0.")

        x0 = int(math.floor(centerX - roiWidth/2))
        x1 = int(math.ceil(centerX + roiWidth/2))
        y0 = int(math.floor(centerY - roiHeight/2))
        y1 = int(math.ceil(centerY + roiHeight/2))

        if x1<0 or y1<0:
            raise Exception("Right or lower boundary for ROI (x1 or y1) cannot be below zero.")

        if roiWidth>self.getWidth() or roiHeight>self.getHeight():
            raise Exception("Size of the ROI is bigger than the image size. ROI: " + str(roiWidth) + " x " + str(roiHeight) + ". Image: " + str(self.getWidth()) + " x " + str(self.getHeight()))   
        if x0 < 0:
            x1 += abs(x0)
            x0 = 0

        if y0 < 0:
            y1 += abs(y0)
            y0 = 0

        if x1 >= self.getWidth():
            x1 = self.getWidth()
            x0 = x1 - roiWidth

        if y1 >= self.getHeight():
            y1 = self.getHeight()
            y0 = y1 - roiHeight

        # These should match roiWidth and roiHeight...
        roiDimX = x1 - x0
        roiDimY = y1 - y0

        self.crop(x0, y0, x1, y1)
        return x0, x1, y0, y1

    def bin(self, binSizeX, binSizeY, operation="mean"):
        """ Decrease image size by merging pixels using specified operation.
            Valid operations: mean, max, min, sum. """

        if binSizeX is None:
            binSizeX = 1

        if binSizeY is None:
            binSizeY = 1

        if (binSizeX > 1) or (binSizeY > 1):
            # Picture dimensions must be integer multiple of binning factor. If not, crop:
            overhangX = math.fmod(int(self.getWidth()), binSizeX)
            overhangY = math.fmod(int(self.getHeight()), binSizeY)
            if (overhangX > 0) or (overhangY > 0):
                #log("Cropping before binning because of nonzero overhang: (" + str(overhangX) + ", " + str(overhangY) + ")")
                self.crop(0, 0, self.getWidth()-int(overhangX), self.getHeight()-int(overhangY))

            newWidth  = self.width // binSizeX
            newHeight = self.height // binSizeY

            # Shift pixel values that need to be binned together into additional axes:
            binshape = (newHeight, binSizeY, newWidth, binSizeX)
            self.px = self.px.reshape(binshape)
            
            # Perform binning operation along binning axes (axis #3 and #1).
            # These axes will be collapsed to contain only the result
            # of the binning operation.
            if operation == "mean":
                self.px = self.px.mean(axis=(3, 1))
            elif operation == "sum":
                self.px = self.px.sum(axis=(3, 1))
            elif operation == "max":
                self.px = self.px.max(axis=(3, 1))
            elif operation == "min":
                self.px = self.px.min(axis=(3, 1))
            elif operation is None:
                raise Exception("No binning operation specified.")
            else:
                raise Exception("Invalid binning operation: {}.".format(operation))

            self.setWidth(newWidth)
            self.setHeight(newHeight)

            # Resolution assumes isotropic pixels...
            self.resolution *= binSizeX

    def addImage(self, other):
        """ Add pixel values from another image to this image. """
        if self.dimensionsMatch(other):
            self.px = self.px + other.getPixelMap()

    def subtractImage(self, other):
        """ Subtract pixel values of another image from this image. """
        if self.dimensionsMatch(other):
            self.px = self.px - other.getPixelMap()

    def multiplyImage(self, other):
        """ Multiply pixel values from another image to this image. """
        if self.dimensionsMatch(other):
            self.px = self.px * other.getPixelMap()

    def divideImage(self, other):
        """ Multiply pixel values by another image. """
        if self.dimensionsMatch(other):
            self.px = self.px / other.getPixelMap()

    def square(self):
        self.px *= self.px

    def sqrt(self):
        self.px = numpy.sqrt(self.px)

    def add(self, value):
        self.px += value

    def subtract(self, value):
        self.px -= value

    def multiply(self, value):
        self.px *= value

    def divide(self, value):
        """ Divide all pixels values by given scalar value. """
        self.px = self.px / float(value)

    def invert(self, min=0, maximum=65535):
        self.px = maximum - self.px

    def renormalize(self, newMin=0, newMax=1, currentMin=None, currentMax=None, ROI=None):
        """Renormalization of grey values from (currentMin, Max) to (newMin, Max) """

        # Take full image if no ROI is given
        if ROI==None:
            ROI = ImageROI(0, 0, self.getWidth(), self.getHeight())

        slc = self.px[ROI.y0:ROI.y1, ROI.x0:ROI.x1]

        if currentMin is None:
            currentMin = slc.min()

        if currentMax is None:
            currentMax = slc.max()

        if(currentMax != currentMin):
            slc = (slc-currentMin)*(newMax-newMin)/(currentMax-currentMin)+newMin
            self.px[ROI.y0:ROI.y1, ROI.x0:ROI.x1] = slc
        else:
            slc = slc*0
            self.px[ROI.y0:ROI.y1, ROI.x0:ROI.x1] = slc
            #raise Exception("Division by zero upon renormalization: currentMax=currentMin={}".format(currentMax))

    def map_lookup(self, gv, gv_from, gv_to):
        """ Return new grey value for given grey value 'gv'. Helper function for self.map()."""

        if gv in gv_from:
            # Given grey value is defined in 'from' list:
            return gv_to[numpy.where(gv_from==gv)]
        else:
            # Linear interpolation:
            a = 0  # left index of interpolation region
            if len(gv_from) > 2:
                for i in range(len(gv_from)-2):
                    if gv_from[i+1] > gv:
                        break

                    a += 1

            b = a + 1  # right index of interpolation region

            xa = gv_from[a]
            xb = gv_from[b]
            ya = gv_to[a]
            yb = gv_to[b] 

            # Slope of linear function:
            m = (yb-ya) / (xb-xa)

            # y axis intersection point ("offset"):
            n = yb - m*xb

            # newly assigned grey value:
            return (m*gv + n)


    def map(self, gv_from, gv_to, bins=1000):
        """ Applies a lookup table (LUT map) to convert image grey values
            according to given assignment tables (two numpy lists).

            gv_from: numpy array of given grey values (in current image)
            gv_to:   numpy array of assigned grey values (for converted image)

            Linear interpolation will take place for gaps in lookup table.
        """

        if len(gv_from) == len(gv_to):
            if len(gv_from) > 1:
                gvMin = self.min()
                gvMax = self.max()

                # Left position of each bin:
                positions, gvStepsize = numpy.linspace(start=gvMin, stop=gvMax, num=bins+1, endpoint=True, dtype=numpy.float64, retstep=True)

                # New grey value for each left position:
                mappingFunction = numpy.vectorize(pyfunc=self.map_lookup, excluded={1, 2})
                newGV = mappingFunction(positions, gv_from, gv_to)

                # Differences in newGV:
                deltaGV = numpy.diff(newGV, n=1)


                # Prepare parameters m (slope) and n (offset) for linear
                # interpolation functions of each bin:
                slopes  = numpy.zeros(bins, dtype=numpy.float64)
                offsets = numpy.zeros(bins, dtype=numpy.float64)

                slopes  = deltaGV / gvStepsize

                #print("newGV:     {}".format(numpy.shape(newGV)))
                #print("slopes:    {}".format(numpy.shape(slopes)))
                #print("positions: {}".format(numpy.shape(positions)))

                offsets = newGV[1:] - slopes*positions[1:]

                inverse_stepsize = 1.0 / gvStepsize

                maxIndices = numpy.full(shape=numpy.shape(self.px), fill_value=bins-1, dtype=numpy.uint32)
                bin_indices = numpy.minimum(maxIndices, numpy.floor((self.px - gvMin) * inverse_stepsize).astype(numpy.uint32))

                m_px = slopes[bin_indices]
                n_px = offsets[bin_indices]

                self.px = m_px*self.px + n_px
            else:
                raise Exception("image.map(): At least two mappings are required in the grey value assignment lists.")
        else:
            raise Exception("image.map(): gv_from must have same length as gv_to.")

    def stats(self, ROI=None):
        """ Image or ROI statistics. Mean, Standard Deviation """

        # Take full image if no ROI is given
        if ROI==None:
            ROI = ImageROI(0, 0, self.getWidth(), self.getHeight())

        slc = self.px[ROI.y0:ROI.y1, ROI.x0:ROI.x1]

        mean  = numpy.mean(slc)
        sigma = numpy.std(slc)
        snr   = 0
        if sigma > 0:
            snr = mean / sigma

        return {"mean": mean, "stddev": sigma, "snr": snr, "width": ROI.width(), "height": ROI.height(), "area": ROI.area()}

    def noise(self, sigma):
        """ Add noise to image.

            Gaussian noise:
            sigma: standard deviation (scalar or array that matches image size)
        """

        rng = default_rng()
        self.px += rng.normal(loc=0, scale=sigma, size=numpy.shape(self.px))

    def smooth_gaussian(self, sigma):
        self.px = ndimage.gaussian_filter(input=self.px, sigma=sigma, order=0, )

    def applyMedian(self, kernelSize=1):
        if kernelSize > 1:
            self.px = ndimage.median_filter(self.px, int(kernelSize))

    def applyThreshold(self, threshold, lower=0, upper=65535):
        self.px = numpy.where(self.px > threshold, upper, lower).astype(self.getInternalDataType())

    def renormalizeToMeanAndStdDev(self, mean, stdDev, ROI=None):
        """ Renormalize grey values such that mean=30000, (mean-stdDev)=0, (mean+stdDev)=60000 """

        # Take full image if no ROI is given
        if ROI==None:
            ROI = ImageROI(0, 0, self.getWidth(), self.getHeight())

        self.px[ROI.y0:ROI.y1, ROI.x0:ROI.x1] = ((self.px[ROI.y0:ROI.y1, ROI.x0:ROI.x1] - mean)/stdDev)*30000 + 30000

    def edges_sobel(self):
        # Sobel edge detection:
        edgesX = ndimage.sobel(self.px, axis=0, mode='nearest')
        edgesY = ndimage.sobel(self.px, axis=1, mode='nearest')
        return numpy.sqrt(edgesX**2 + edgesY**2)

    def edges_canny(self):
        # The 'feature' package from scikit-image,
        # only needed for Canny edge detection, when used instead of Sobel.
        from skimage.feature import canny   # Canny edge detection

        # Canny edge detection. Needs 'scikit-image' package.  from skimage import feature
        return canny(self.px)

    def filter_edges(self, mode='sobel'):
        if(mode == 'sobel'):
            self.px = self.edges_sobel()
        elif(mode == 'canny'):
            self.px = self.edges_canny()
        else:
            raise Exception("Valid edge detection modes: 'sobel'")
        
        # Rescale:
        self.px = self.px.astype(self.getInternalDataType())
        #self.thresholding(0)    # black=0, white=65535

    def cleanPatches(self, min_patch_area=None, max_patch_area=None, remove_border_patches=False, aspect_ratio_tolerance=None):
        iterationStructure = ndimage.generate_binary_structure(rank=2, connectivity=2)  # apply to rank=2D array, only nearest neihbours (connectivity=1) or next nearest neighbours as well (connectivity=2)

        labelField, nPatches = ndimage.label(self.px, iterationStructure)
        nCleaned   = 0
        nRemaining = 0
        patchGeometry = []

        if nPatches == 0:
            log("Found no structures")
        else:
            self.erase()

            areaMin = 0
            if(min_patch_area != None):
                areaMin = min_patch_area
            
            areaMax = self.getWidth() * self.getHeight()
            if(max_patch_area != None):
                areaMax = max_patch_area

            areaMin = areaMin / (self.getResolution()**2)
            areaMax = areaMax / (self.getResolution()**2)

            for i in range(1, nPatches+1):
                patchCoordinates = numpy.nonzero(labelField==i)

                # Check patch size:
                nPatchPixels = len(patchCoordinates[0])
                if nPatchPixels < areaMin or nPatchPixels > areaMax:  # Black out areas that are too small or too big for a circle
                    nCleaned += 1
                    continue
                
                coordinatesX = patchCoordinates[1]
                coordinatesY = patchCoordinates[0]

                left  = numpy.amin(coordinatesX)
                right = numpy.amax(coordinatesX)
                top   = numpy.amin(coordinatesY)
                bottom= numpy.amax(coordinatesY)

                if remove_border_patches:   
                    if((left==0) or (top==0) or (right==self.getWidth()-1) or (bottom==self.getHeight()-1)):
                        nCleaned += 1
                        continue

                # An ideal circle should have an aspect ratio of 1:
                if aspect_ratio_tolerance != None:
                    aspectRatio = 0
                    if(top != bottom):
                        aspectRatio = abs(right-left) / abs(bottom-top)

                    if abs(1-aspectRatio) > aspect_ratio_tolerance:  # This is not a circle
                        nCleaned += 1
                        log("Aspect ratio {ar:.3f} doesn't meet aspect ratio tolerance |1-AR|={tolerance:.3f}".format(ar=aspectRatio, tolerance=aspect_ratio_tolerance))
                        continue

                # Add patch center as its coordinate:
                patchGeometry.append(((right+left)/2.0, (bottom+top)/2.0, right-left, bottom-top))

                self.px[patchCoordinates] = 1
                nRemaining += 1

        return nPatches, nCleaned, nRemaining, patchGeometry

    def fitCircle(self):
        # Linear least squares method by:
        # I. D. Coope,
        # Circle Fitting by Linear and Nonlinear Least Squares,
        # Journal of Optimization Theory and Applications, 1993, Volume 76, Issue 2, pp 381-388
        # https://doi.org/10.1007/BF00939613

        coordinates = numpy.nonzero(self.px)
        circlePixelsX = coordinates[1]
        circlePixelsY = coordinates[0]
        nPoints = len(circlePixelsX)
        circlePixels1 = numpy.ones(nPoints)

        # Create the matrix B for the system of linear equations:
        matrixB = numpy.array((circlePixelsX, circlePixelsY, circlePixels1))
        matrixB = matrixB.transpose()

        # linear equation to optimize:
        # matrix B * result = vector d
        d = []
        for i in range(nPoints):
            d.append(circlePixelsX[i]**2 + circlePixelsY[i]**2)

        vectorD = numpy.array(d)

        results, residuals, rank, s = numpy.linalg.lstsq(matrixB, vectorD, rcond=None)

        centerX = (results[0] / 2.0)
        centerY = (results[1] / 2.0)
        radius  = math.sqrt(results[2] + centerX**2 + centerY**2)

        # Calculate deviation statistics:
        differenceSum = 0
        minDifference = 99999
        maxDifference = 0
        for i in range(nPoints):
            diff = abs(radius  -  math.sqrt((centerX - circlePixelsX[i])**2 + (centerY - circlePixelsY[i])**2))
            differenceSum += diff

            if minDifference > diff:
                minDifference = diff

            if maxDifference < diff:
                maxDifference = diff

        meanDifference = differenceSum / nPoints

        return centerX, centerY, radius, meanDifference, minDifference, maxDifference

    def intensityFunction2D(self, x, I0, mu, R, x0):   # Lambert-Beer-Law for ball intensity, to fit.
        radicand = numpy.power(R,2) - numpy.power((x-x0),2)
        
        # Avoid root of negative numbers
        radicand[radicand < 0] = 0   

        # Huge radicands lead to exp()->0, therefore avoid huge exponentiation:
        radicand[radicand > (1400*1400)] = (1400*1400)

        result = I0*numpy.exp(-2.0*mu*numpy.sqrt(radicand))

        return result

    def intensityFunction3D(self, coord, I0, mu, R, x0, y0):   # Lambert-Beer-Law for ball intensity, to fit.
        if len(coord) == 2:
            (x, y) = coord

            radicand = numpy.power(R,2) - numpy.power((x-x0),2) - numpy.power((y-y0),2)
            
            # Avoid root of negative numbers
            radicand[radicand < 0] = 0   

            # Huge radicands lead to exp()->0, therefore avoid huge exponentiation:
            radicand[radicand > (1400*1400)] = (1400*1400)

            result = I0 * numpy.exp(-2.0*mu*numpy.sqrt(radicand))
            
            return result
        else:
            raise Exception("3D Intensity fit function expects a tuple (x,y) for coordinates.")

    def fitIntensityProfile(self, axis="x", initI0=None, initMu=0.003, initR=250, initX0=None, avgLines=5):
        yData = 0
        xdata = 0
        if initI0 is None:
            initI0 = self.max()   # Hoping that a median has been applied before.

        if axis == "x":
            if initX0 is None:
                initX0 = self.getWidth() / 2

            startLine = int((self.getHeight() / 2) - math.floor(avgLines/2))
            stopLine  = int((self.getHeight() / 2) + math.floor(avgLines/2))

            # Accumulate intensity profile along 'avgLines' lines around the center line:
            yData = numpy.zeros(self.getWidth(), dtype=self.getInternalDataType())
            for l in range(startLine, stopLine+1):
                yData += self.px[l,:]

            xData = numpy.linspace(0, self.getWidth()-1, self.getWidth())

        elif axis == "y":
            if initX0 is None:
                initX0 = self.getHeight() / 2

            startLine = int((self.getWidth() / 2) - math.floor(avgLines/2))
            stopLine  = int((self.getWidth() / 2) + math.floor(avgLines/2))

            # Accumulate intensity profile along 'avgLines' lines around the center line:
            yData = numpy.zeros(self.getHeight(), dtype=self.getInternalDataType())
            for l in range(startLine, stopLine+1):
                yData += self.px[:,l]

            xData = numpy.linspace(0, self.getHeight()-1, self.getHeight())

        else:
            raise Exception("projectionImage::fitIntensityProfile() needs profile direction to be 'x' or 'y'.")

        yData = yData / int(avgLines)   # average intensity profile
        firstGuess = (initI0, initMu, initR, initX0)

        try:
            optimalParameters, covariances = optimize.curve_fit(self.intensityFunction2D, xData, yData, p0=firstGuess)
        except Exception:
            optimalParameters = (None, None, None, None)


        fittedI0 = optimalParameters[0]
        fittedMu = optimalParameters[1]
        fittedR  = optimalParameters[2]
        fittedX0 = optimalParameters[3]

        return fittedI0, fittedMu, fittedR, fittedX0

class ImageStack:
    """ Specify an image stack from a single file (RAW chunk) or
        a collection of single 2D RAW or TIFF files. """

    def __init__(self, filePattern=None, width=None, height=None, dataType=None, byteOrder=None, rawFileHeaderSize=0, rawImageHeaderSize=0, slices=None, startNumber=0, flipByteOrder=False):
        self.files = ImageFile(filePattern, dataType, byteOrder, flipByteOrder)

        # Has this stack already been built?
        self.built = False

        self.width       = width
        self.height      = height
        self.nSlices     = slices   # number of slices in stack
        self.startNumber = startNumber

        # A RAW chunk can contain an overall file header, and
        # each image in the stack can contain an image header.
        self.rawFileHeaderSize = rawFileHeaderSize
        self.rawImageHeaderSize = rawImageHeaderSize

        self._isVolumeChunk = False    # Is this a volume chunk or is a file list provided?

        self.fileList = []
        self.fileNumbers = []   # store original stack number in file name

    def addStack(self, other):
        if (self.width == other.width) and (self.height == other.height):
            self.nSlices += other.nSlices
            self.fileList.extend(other.fileList)
            self.fileNumbers.extend(other.fileNumbers)
        else:
            raise Exception("Error adding stack: image dimensions don't match.")

    def isVolumeChunk(self):
        return self._isVolumeChunk

    def setVolumeChunk(self, isVolumeChunk):
        self._isVolumeChunk = isVolumeChunk

    def getFileByteOrder(self):
        return self.files.getByteOrder()

    def setFileByteOrder(self, byteOrder):
        self.files.setByteOrder(byteOrder)

    def getFileDataType(self):
        return self.files.getDataType()

    def setFileDataType(self, dataType):
        self.files.setDataType(dataType)

    def doFlipByteOrder(self):
        return self.files.doFlipByteOrder()

    def setFlipByteOrder(self, flipByteOrder):
        self.files.setFlipByteOrder(flipByteOrder)

    def fileStackInfo(self, filenameString):
        """ Split file pattern into lead & trail text, number of expected digits. """
        if '%' in filenameString:
            # A % sign in the provided file pattern indicates an image stack: e.g. %04d
            percentagePosition = filenameString.find("%")

            numberStart = percentagePosition + 1
            numberStop  = filenameString.find("d", percentagePosition)

            leadText  = ""
            if(percentagePosition > 0):
                leadText = filenameString[:percentagePosition]

            trailText = ""
            if((numberStop+1) < len(filenameString)):
                trailText = filenameString[(numberStop+1):]

            if(numberStop > numberStart):
                numberString = filenameString[numberStart:numberStop]
                if(numberString.isdigit()):
                    nDigitsExpected = int(numberString)
                    return leadText, trailText, nDigitsExpected
                else:
                    raise Exception("Image stack pattern is wrong. The wildcard for sequential digits in a filename must be %, followed by number of digits, followed by d, e.g. %04d")
            else:
                raise Exception("Image stack pattern is wrong. The wildcard for sequential digits in a filename must be %, followed by number of digits, followed by d, e.g. %04d")

        return filenameString, "", 0

    def buildStack(self):
        """ Build list of files that match given file name pattern. """
        self.fileList = []
        self.fileNumbers = []

        # Treat projection files
        inFilePattern = self.files.getFilename()
        inputFolder  = os.path.dirname(inFilePattern)
        projBasename = os.path.basename(inFilePattern)

        if inputFolder == "" or inputFolder is None:
            inputFolder = "."

        # Check if an image stack is provided:
        if('%' not in inFilePattern):
            self.fileList.append(inFilePattern)

            if(isTIFF(inFilePattern)):  # treat as single TIFF projection            
                self._isVolumeChunk = False
                testImage = Image(inFilePattern)
                testImage.read()
                self.width    = testImage.getWidth()
                self.height   = testImage.getHeight()
                self.nSlices  = 1
                self.files.setDataType(testImage.inputFile.getDataType())
            else:  # treat as raw chunk
                if (self.width != None) and (self.height != None):
                    if (self.files.getDataType() != None):
                        if os.path.isfile(inFilePattern):
                            self._isVolumeChunk = True

                            if (self.nSlices is None):
                                # Determine number of slices.
                                fileSizeInBytes = os.path.getsize(inFilePattern)
                                dataSizeInBytes = fileSizeInBytes - self.rawFileHeaderSize
                                bytesPerImage = self.rawImageHeaderSize + self.width * self.height * self.files.getDataType().itemsize

                                if (dataSizeInBytes >= bytesPerImage):
                                    if (dataSizeInBytes % bytesPerImage) == 0:
                                        self.nSlices = int(dataSizeInBytes / bytesPerImage)
                                        log("{} slices found in raw chunk.".format(self.nSlices))
                                    else:
                                        raise Exception("The raw chunk data size ({} bytes, without general file header) is not divisible by the calculated size of a single image ({} bytes, including image header). Therefore, the number of slices cannot be determined. {}".format(dataSizeInBytes, bytesPerImage, inFilePattern))
                                else:
                                    raise Exception("The raw chunk data size ({} bytes, without general file header) is smaller than the calculated size of a single image ({} bytes, including image header). {}".format(dataSizeInBytes, bytesPerImage, inFilePattern))
                        else:
                            raise Exception("File not found: {}".format(inFilePattern))
                    else:
                        raise Exception("Please provide the data type of the raw chunk.")
                else:
                    raise Exception("Please provide width and height (in pixels) of the raw chunk.")
        else:
            # A % sign in the provided file pattern indicates an image stack: e.g. %04d
            leadText, trailText, nDigitsExpected = self.fileStackInfo(projBasename)

            # Get list of files in input folder:
            fileList = os.listdir(inputFolder)
            fileList.sort()

            nImported = 0

            for f in fileList:
                file = inputFolder + "/" + f
                if os.path.isfile(file):
                    # Check if filename matches pattern:
                    if(f.startswith(leadText) and f.endswith(trailText)):
                        digitText = f[len(leadText):-len(trailText)]
                        if digitText.isdigit(): # and len(digitText)==nDigitsExpected:
                            # Pattern matches.
                            n = int(digitText)
                            if n >= self.startNumber:
                                self.fileList.append(file)
                                self.fileNumbers.append(n)

                                nImported += 1
                                if nImported == self.nSlices:
                                    break
                        else:
                            continue
                    else:
                        continue

            self.nSlices = len(self.fileList)

            if self.nSlices > 0:
                if isTIFF(self.fileList[0]):
                    testImage = Image(self.fileList[0])
                    testImage.read()
                    self.width    = testImage.getWidth()
                    self.height   = testImage.getHeight()
                    self.files.setDataType(testImage.inputFile.getDataType())

        self.built = True
                

    def getFilename(self, index=None):
        if index != None:
            if self._isVolumeChunk:
                if len(self.fileList) > 0:
                    return self.fileList[0]
                else:
                    return None
            else:
                if len(self.fileList) > index:
                    return self.fileList[index]
                else:
                    return None
        else:
            return self.files.getFilename()

    def getFileBasename(self, index=None):
        if index != None:
            if self._isVolumeChunk:
                if len(self.fileList) > 0:
                    return os.path.basename(self.fileList[0])
                else:
                    return None
            else:
                if len(self.fileList) > index:
                    return os.path.basename(self.fileList[index])
                else:
                    return None
        else:
            return self.files.getFileBasename()

    def setFilename(self, filename):
        self.files.setFilename(filename)

    def getImage(self, index, outputFile=None):
        """ Read and return image at position 'index' within the stack. """
        if index >= 0:
            if not self._isVolumeChunk:  # read single image file from stack:
                if len(self.fileList) > index:
                    filename = self.fileList[index]
                    file = ImageFile(filename=filename, dataType=self.getFileDataType(), byteOrder=self.getFileByteOrder(), flipByteOrder=self.doFlipByteOrder())

                    img = Image(file, outputFile)
                    if isTIFF(filename):
                        img.read()
                    else:
                        img.readRAW(self.width, self.height, 0, self.getFileDataType(), self.getFileByteOrder(), self.rawFileHeaderSize, self.rawImageHeaderSize)
                    return img
                else:
                    raise Exception("The requested slice nr. {} is out of bounds, because only {} image files were found.".format(index, len(self.fileList)))
            else:  # read slice from volume chunk, obeying start number
                if len(self.fileList) > 0:
                    file = self.fileList[0]
                    img = Image(file, outputFile)
                    chunkIndex = index + self.startNumber
                    if isTIFF(file):
                        raise Exception("Cannot treat 3D TIFFs.")
                    else:
                        img.readRAW(self.width, self.height, chunkIndex, self.getFileDataType(), self.getFileByteOrder(), self.rawFileHeaderSize, self.rawImageHeaderSize)
                        return img
                else:
                    raise Exception("No image file specified to be loaded.")
        else:
            raise Exception("Negative slice numbers do not exists. {} requested.".format(index))

    def getMeanImage(self, outputFile=None):
        """ Calculate the mean of all image files. """
        if self.nSlices > 0:
            if self.nSlices > 1:
                sumImg = self.getImage(0, outputFile)
                for i in range(1, self.nSlices):
                    print("\rMean Image: summing up {i}/{n}".format(i=(i+1), n=self.nSlices), end='')
                    sumImg.addImage(self.getImage(i, outputFile))
                    

                print("")

                sumImg.divide(self.nSlices)
                return sumImg
            else:
                return self.getImage(0, outputFile)
        else:
            return None

    def getStdDevImage(self, meanImg=None, outputFile=None):
        """ Calculate the pixel-wise RMS of the image files. """
        if self.nSlices > 0:
            if self.nSlices > 1:
                if meanImg is None:
                    meanImg = self.getMeanImage(outputFile)

                sumImg = Image()
                sumImg.shapeLike(otherImg=meanImg)

                for i in range(0, self.nSlices):
                    print("\rRMSD Image: component {i}/{n}".format(i=i+1, n=self.nSlices), end='')
                    sqDiffImg = self.getImage(i, outputFile)
                    sqDiffImg.subtractImage(meanImg)
                    sqDiffImg.square()

                    sumImg.addImage(sqDiffImg)

                sumImg.divide(self.nSlices)
                sumImg.sqrt()

                print("")

                return sumImg
            else:
                return self.getImage(0, outputFile)
        else:
            return None