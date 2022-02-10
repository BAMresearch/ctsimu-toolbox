# -*- coding: UTF-8 -*-

# Version 0.1 (2020-09-29)

import os      # to check for files
import sys     # for native byteorder
import struct  # for C-style byte reading
import math    # for floor
import copy    # for deepcopy
import numpy   # for image data arrays

# Define TIFF field types:
TIFF_BYTE      =  1  #  8 bit unsigned
TIFF_ASCII     =  2  #  8 bit character
TIFF_SHORT     =  3  # 16 bit unsigned
TIFF_LONG      =  4  # 32 bit unsigned
TIFF_RATIONAL  =  5  # 2x32 bit unsigned (numerator, denominator)
TIFF_SBYTE     =  6  #  8 bit signed
TIFF_UNDEFINED =  7  #  8 bit of unknown data
TIFF_SSHORT    =  8  # 16 bit signed
TIFF_SLONG     =  9  # 32 bit signed
TIFF_SRATIONAL = 10  # 2x32 bit signed (numerator, denominator)
TIFF_FLOAT     = 11  #  4 byte single precision IEEE float
TIFF_DOUBLE    = 12  #  8 byte double precision IEEE float

TIFF_BYTES_PER_TYPE = [1, 1, 2, 4, 8, 1, 1, 2, 4, 8, 4, 8]
TIFF_STRUCT_CHAR_FOR_TYPE = ["B", "c", "H", "L", "LL", "b", "B", "h", "l", "ll", "f", "d"]

TIFF_INTTYPES = [TIFF_BYTE, TIFF_SHORT, TIFF_LONG, TIFF_SBYTE, TIFF_SSHORT, TIFF_SLONG]

# TIFF values
# Photometric Interpretation
TIFF_WHITE_IS_ZERO = 0
TIFF_BLACK_IS_ZERO = 1
TIFF_RGB           = 2

# Compression
TIFF_NO_COMPRESSION       = 1
TIFF_CCITT_COMPRESSION    = 2
TIFF_PACKBITS_COMPRESSION = 32773
TIFF_LZW_COMPRESSION      = 5

# Resolution unit
TIFF_RES_NONE       = 1
TIFF_RES_INCH       = 2
TIFF_RES_CENTIMETER = 3

# Sample format
TIFF_SAMPLEFORMAT_UINT   = 1
TIFF_SAMPLEFORMAT_INT    = 2
TIFF_SAMPLEFORMAT_IEEEFP = 3
TIFF_SAMPLEFORMAT_VOID   = 4

# Planar configuration
TIFF_CHUNKY = 1
TIFF_PLANAR = 2

# Fill order
TIFF_MSB2LSB = 1
TIFF_LSB2MSB = 2

# Predictor
TIFF_NO_PREDICTOR = 1
TIFF_HORIZONTAL_DIFFERENCING = 2

TAGNAMES = {
    "254": "New subfile type",
    "256": "Cols",
    "257": "Rows",
    "258": "Bits per sample",
    "259": "Compression",
    "262": "Photometric interpretation",
    "266": "Fill order",
    "270": "Image description",
    "273": "Strip offsets",
    "274": "Orientation",
    "277": "Samples per pixel",
    "278": "Rows per strip",
    "279": "Strip byte counts",
    "282": "Resolution X",
    "283": "Resolution Y",
    "284": "Planar configuration",
    "296": "Resolution unit",
    "317": "Predictor",
    "320": "Color map",
    "339": "Sample format"
}

TIFFY_LOGLEVEL_ERROR   = 0
TIFFY_LOGLEVEL_WARNING = 10
TIFFY_LOGLEVEL_INFO    = 20
TIFFY_LOGLEVEL_DEBUG   = 30

tiffy_currentLogLevel = TIFFY_LOGLEVEL_ERROR

def tiffyLog(level, message):
    if level < tiffy_currentLogLevel:
        print(message)


def getByteOrder(numpyArray):
    byteOrder = numpyArray.dtype.byteorder
    if byteOrder == "=":  # native byte order
        if sys.byteorder == "big":
            byteOrder = ">"
        else:
            byteOrder = "<"

    return byteOrder

def tagName(intTag):
    strTag = "{}".format(intTag)
    if strTag in TAGNAMES:
        return TAGNAMES[strTag]
    else:
        return "Unknown"

def bytesPerTIFFtype(tp):
    if tp > 0 and tp <= len(TIFF_BYTES_PER_TYPE):
        return TIFF_BYTES_PER_TYPE[tp-1]

    raise Exception("Unknown TIFF data type: {}".format(tp))

def TIFFtypeStructCharacter(tp):
    if tp > 0 and tp <= len(TIFF_BYTES_PER_TYPE):
        return TIFF_STRUCT_CHAR_FOR_TYPE[tp-1]

    raise Exception("Unknown TIFF data type: {}".format(tp))

class bits:
    def __init__(self):
        self.data = bytearray()

    def set(self, n):
        """ Takes a number n and stores it in self.data in bitwise manner. """
        self.data = bytearray()

        pos = 0
        while pos>0:
            self.setBit(pos, n&1)
            n = n >> 1
            pos += 1

    def __str__(self):
        """ Print in binary representation. """
        maxPos = len(self.data)*8
        binString = ""

        pos = 0
        nChunks = 0
        nBytes = 0
        while pos <= maxPos:
            if pos%9 == 0:
                binString = "." + binString
                nChunks += 1
            if pos%8 == 0:
                binString = "|" + binString
                nBytes += 1

            binString = "{}".format(self.getBit(pos)) + binString
            pos += 1

        binString += "\n{} chunks in {} bytes.".format(nChunks, nBytes)
        return binString

    def setBytes(self, b):
        """ Takes a bytes object b and stores it in self.data. """
        self.data = bytearray()
        self.data += b

    def reverseBitsInBytes(self):
        for i in range(len(self.data)):
            newByte = 0
            for b in range(8):
                newByte = newByte << 1
                newByte += ((self.data[i] >> b) & 1)
                
            self.data[i] = newByte

    """
    def reverseBits(self, start, stop):
        newBits = bits()
        pos = 0
        for i in range(stopBit-1, startBit-1, -1):
            newBits.setBit(pos, self.getBit(i))
            pos += 1

        for i in range(startBit, stopBit):
            self.setBit(i, self.getBit(i-startBit))
    """

    def byteAndBit(self, pos):
        """ From a given bit position, returns the byte index and the bit index within this byte. """
        byteIdx = int(math.floor(pos/8))
        inByte = pos%8

        return byteIdx, inByte

    def getInt(self, startBit, stopBit):
        s = 0
        if startBit < stopBit:
            for i in range(stopBit-1, startBit-1, -1):
                s = s << 1
                s += self.getBit(i)

        return int(s)

    def getIntMSBtoLSB(self, startBit, stopBit):
        s = 0
        for i in range(startBit, stopBit):  # MSB to LSB
            s = s << 1
            s += self.getBit(i)

        return s

    def getIntMSBtoLSB_inBytes(self, startByte, startBit, nBits):
        s = 0

        # Maximum of three masks required (for 9..12 bits)
        mask1 = mask << startBit
        mask2 = 255 << (nBits-7)
        mask3 = 255

        """
        i = 0
        #print("StartByte: {}, StartBit: {}, nBits: {}".format(startByte, startBit, nBits))
        while i<nBits:  # MSB to LSB
            s = s << 1  # Shift left by 1
            s += ((self.data[startByte] >> (7-startBit)) & 1)

            startBit += 1
            if startBit > 7:
                startBit = 0
                startByte += 1

            i += 1
        """

        startBit += nBits
        if startBit > 7:
            startBit = 0
            startByte += 1

        return s, startByte+1, startBit+nBits

    def getIntMSBtoLSB_faster(self, startByte, rightZeros, mask0, mask1, mask2):
        s = ((self.data[startByte]&mask0)<<16) + ((self.data[startByte+1]&mask1)<<8) + (self.data[startByte+2]&mask2)

        s = s >> rightZeros
        #newByte = 0
        #for b in range(9):
        #    newByte = newByte << 1
        #    newByte += ((s >> b) & 1)

        #return newByte

        return s

    def getBit(self, pos):
        byteIdx, inByte = self.byteAndBit(pos)

        if(byteIdx >= len(self.data)):
            return 0

        return ((self.data[byteIdx] >> inByte) & 1)

    def getBitInByte(self, byteIdx, bitPos):
        return ((self.data[byteIdx] >> bitPos) & 1)

    def setBit(self, pos, value=1):
        byteIdx, inByte = self.byteAndBit(pos)

        while(byteIdx > len(self.data)):
            self.data += bytes(1)

        mask = 1 << inByte
        if value == 1:   # set bit
            self.data[byteIdx] = self.data[byteIdx] | mask
        else:  # clear bit
            self.data[byteIdx] = self.data[byteIdx] & ~mask

    def setBitInByte(self, byteIdx, bitPos, value=1):
        while(byteIdx > len(self.data)):
            self.data += bytes(1)

        mask = 1 << bitPos
        if value == 1:   # set bit
            self.data[byteIdx] = self.data[byteIdx] | mask
        else:  # clear bit
            self.data[byteIdx] = self.data[byteIdx] & ~mask       

class lzwStringTable:
    def __init__(self):
        self.byteStrings = []
        self.init()

    def init(self):
        self.byteStrings = []

        # Create all bytes:
        for i in range(256):
            self.byteStrings.append(bytearray(struct.pack("B", i)))

        self.byteStrings.append(bytearray())  # ClearCode: 256
        self.byteStrings.append(bytearray())  # EndOfInformation code: 257

    def currentCodeBitWidth(self):
        """
        l = len(self.byteStrings) + 1
        bitWidth = 0
        while l != 0:
            l = l >> 1
            bitWidth += 1

        return bitWidth
        """

        if len(self.byteStrings) < 255:
            return 8
        elif len(self.byteStrings) < 511:
            return 9
        elif len(self.byteStrings) < 1023:
            return 10
        elif len(self.byteStrings) < 2047:
            return 11
        elif len(self.byteStrings) < 4095:
            return 12
        elif len(self.byteStrings) < 8191:
            return 13
        elif len(self.byteStrings) < 16383:
            return 14
        else:
            raise Exception("LZW: Dictionary is too big.")

    def contains(self, code):
        if code < len(self.byteStrings):
            return True

        return False

    def add(self, b):
        self.byteStrings.append(b)

    def isClearCode(self, code):
        if code == 256:
            return True

        return False

    def isEndOfInformation(self, code):
        if code == 257:
            return True

        return False

    def stringFromCode(self, code, codeWidth):
        #if code < 0:
        #    return bytearray()

        return self.byteStrings[code]

        """
        if code < len(self.byteStrings):
            #print("String from code {}: {}".format(code, self.byteStrings[code]))
            return self.byteStrings[code]
        else:
            raise Exception("LZW: Requested code #{} not in current string table. Current string table size: {}. Current bit width: {}".format(code, len(self.byteStrings), codeWidth))
        """

class lzwData:
    def __init__(self):
        self.compressed = None
        self.umcompressed = None
        self.stringTable = lzwStringTable()

        self.currentBitPos = 0

        self.currentByteIdx = 0
        self.currentBitPosInByte  = 0  # in byte

        self.codeWidth = 9  # Starting with 9 bit wide codes

        self.masks = {} #[[[0]*3]*8]*14
        self.byteSkips = {} # [[0]*8]*14
        self.leftAlignOffsets = {}

        for codeWidth in range(9, 14):
            fundMask = 2**codeWidth - 1

            for shift in range(0, 8):   # Bitshift / startBit
                leftAlignOffset = 3*8 - codeWidth - shift

                mask = fundMask << leftAlignOffset
                mask0 = (mask >> 16) & 255  # 00000000 00000000 11111111
                mask1 = (mask >> 8) & 255   # 00000000 11111111 00000000
                mask2 = (mask) & 255        # 11111111 00000000 00000000

                self.masks[codeWidth,shift,0] = mask0
                self.masks[codeWidth,shift,1] = mask1
                self.masks[codeWidth,shift,2] = mask2

                self.byteSkips[codeWidth,shift] = int(math.ceil((codeWidth+shift+1)/8))-1
                self.leftAlignOffsets[codeWidth,shift] = leftAlignOffset

    def resetStringtable(self):
        self.stringTable.init()
        self.codeWidth = 9

    def setCompressed(self, compressed):
        self.compressed = bits()
        self.compressed.setBytes(compressed)
        self.compressed.data += bytes(3)   # for integer conversion

    def getNextCode(self):
        #s1 = self.currentBitPos
        #s2 = s1 + self.codeWidth

        #code = self.compressed.getIntMSBtoLSB(startBit=s1, stopBit=s2)
        #code, self.currentByteIdx, self.currentBitPosInByte = self.compressed.getIntMSBtoLSB_inBytes(startByte=self.currentByteIdx, startBit=self.currentBitPosInByte, nBits=self.codeWidth)

        startByte = self.currentByteIdx
        mask0=self.masks[self.codeWidth,self.currentBitPosInByte,0]
        mask1=self.masks[self.codeWidth,self.currentBitPosInByte,1]
        mask2=self.masks[self.codeWidth,self.currentBitPosInByte,2]

        code = self.compressed.getIntMSBtoLSB_faster(startByte=startByte, rightZeros=self.leftAlignOffsets[self.codeWidth,self.currentBitPosInByte], mask0=mask0, mask1=mask1, mask2=mask2)

        tiffyLog(TIFFY_LOGLEVEL_DEBUG, "Requesting {}({}). Skip: {} CodeWidth: {} Offset: {}  Masks: {:08b}.{:08b}.{:08b} -> {}".format(startByte, self.currentBitPosInByte, self.byteSkips[self.codeWidth,self.currentBitPosInByte], self.codeWidth, self.currentBitPosInByte, mask0, mask1, mask2, code))

        self.currentByteIdx += self.byteSkips[self.codeWidth,self.currentBitPosInByte]
        self.currentBitPosInByte += self.codeWidth
        self.currentBitPosInByte = self.currentBitPosInByte % 8

        return code

    def decompress(self):
        self.uncompressed = bytearray()

        self.currentBitPos = 0
        oldCode = -1
        code = self.getNextCode()

        self.stringTable.init()  # If compressed data doesn't start with ClearCode

        while not self.stringTable.isEndOfInformation(code):
            if(self.stringTable.isClearCode(code)):
                self.resetStringtable()                

                code = self.getNextCode()
                if(self.stringTable.isEndOfInformation(code)):
                    break

                self.uncompressed += self.stringTable.stringFromCode(code, self.codeWidth)
                oldCode = code
            else:
                if self.stringTable.contains(code):
                    s = self.stringTable.stringFromCode(code, self.codeWidth)
                    self.uncompressed += s
                    self.stringTable.add(self.stringTable.stringFromCode(oldCode, self.codeWidth) + bytes([s[0]]))
                else:
                    s = self.stringTable.stringFromCode(oldCode, self.codeWidth)
                    outString = s + bytes([s[0]])
                    self.uncompressed += outString
                    self.stringTable.add(outString)

                oldCode = code

                self.codeWidth = self.stringTable.currentCodeBitWidth()

            code = self.getNextCode()


class ifdEntry:
    def __init__(self):
        self.tag      = 0  # TIFF Tag ID
        self.type     = 0  # Field Type
        self.count    = 0  # Number of values (NOT bytes)
        self.ifdEntryPos = 0

        self.values = []
        self.extraValuesOffset = 0   # storage offset (in byte) for any extra values that don't fit in the 4 values bytes of this IFD entry

    def set(self, tag, typ, values):
        self.setTagID(tag)
        self.setType(typ)
        self.setValues(values)

    def setTagID(self, tagid):
        self.tag = tagid

    def setType(self, typ):
        self.type = typ

    def setValue(self, value):
        if self.type in TIFF_INTTYPES:
            value = int(value)

        self.values = [value, ]
        self.count = 1

    def setValues(self, values):
        self.values = [0] * len(values)
     
        for i in range(len(self.values)):
            if self.type in TIFF_INTTYPES:
                self.values[i] = int(values[i])
            else:
                self.values[i] = values[i]

        self.count = len(values)

    def nValueBytes(self):
        return self.count * bytesPerTIFFtype(self.type)

    def read(self, offset, f, byteOrder):
        f.seek(offset)
        buff = f.read(8)
        (self.tag, self.type, self.count) = struct.unpack("{endian}HHL".format(endian=byteOrder), buff)

        # Calculate amount of bytes necessary for value(s):
        self.nValueBytes = self.nValueBytes()

        # Read value(s):
        f.seek(offset+8)
        buff = f.read(4)

        if self.nValueBytes > 4:  # The field's actual value bytes must be read from the provided pointer.
            (valueOffset,) = struct.unpack("{endian}L".format(endian=byteOrder), buff)
            f.seek(valueOffset)
            buff = f.read(self.nValueBytes)
      
        # Prepare a struct pattern:
        structPattern = ""
        structCharacter = TIFFtypeStructCharacter(self.type)
        for i in range(self.count):
            structPattern += structCharacter

        # Fill rest of pattern with pad bytes:
        if self.nValueBytes < 4:
            for j in range(4 - self.nValueBytes):
                structPattern += "x"

        tup = struct.unpack("{endian}{pattern}".format(endian=byteOrder, pattern=structPattern), buff)

        if self.type == TIFF_RATIONAL or self.type == TIFF_SRATIONAL:
            for numerator, denominator in zip(tup, tup[1:]):
                val = 0
                if denominator != 0:
                    val = numerator / denominator

                self.values.append(val)
        else:
            for val in tup:
                self.values.append(val)

        if len(self.values) < 5:
            tiffyLog(TIFFY_LOGLEVEL_INFO, " -- Tag {tag} ({tagname}): {n} value(s) {val}".format(tag=self.tag, tagname=tagName(self.tag), n=len(self.values), val=self.values))
        else:
            tiffyLog(TIFFY_LOGLEVEL_INFO, " -- Tag {tag} ({tagname}): {n} value(s)".format(tag=self.tag, tagname=tagName(self.tag), n=len(self.values)))

    def getValue(self):
        if len(self.values) > 0:
            return self.values[0]
        else:
            raise Exception("TIFF field with tag {tag} does not come with any value.".format(tag=self.tag))

    def sizeInBytes(self):
        size = 12
        size += self.sizeOfExtraValues()

        return size

    def sizeOfExtraValues(self):
        size = 0

        # If the values don't fit in the IFD entry, they are stored somewhere else:
        if (bytesPerTIFFtype(self.type)*len(self.values)) > 4:
            size += bytesPerTIFFtype(self.type)*len(self.values)

        return size

    def printOffset(self):
        tiffyLog(TIFFY_LOGLEVEL_INFO, "    IFD entry, tag {t}  {offset}".format(t=self.tag, offset=self.ifdEntryPos))

    def printExtraDataOffset(self):
        if self.nValueBytes() > 4:
            tiffyLog(TIFFY_LOGLEVEL_INFO, "    IFD entry, tag {t}  {offset}  Extra Data".format(t=self.tag, offset=self.extraValuesOffset))

    def prepareDataOffsets(self, offset, extraDataOffset):
        self.ifdEntryPos = offset
        self.extraValuesOffset = extraDataOffset
        return (offset+12), (extraDataOffset + self.sizeOfExtraValues())

    def write(self, f, byteOrder):
        tiffyLog(TIFFY_LOGLEVEL_INFO, "IDF Entry at pos {}".format(f.tell()))
        buff = struct.pack("{endian}HHL".format(endian=byteOrder), self.tag, self.type, self.count)
        f.write(buff)

        if self.nValueBytes() > 4:  # point to value storage offset
            tiffyLog(TIFFY_LOGLEVEL_DEBUG, "Packing for tag {tag} (type {typ}, count {cnt}): Pointing to {offset}".format(tag=self.tag, typ=self.type, cnt=self.count, endian=byteOrder, offset=self.extraValuesOffset))
            buff = struct.pack("{endian}L".format(endian=byteOrder), self.extraValuesOffset)
            f.write(buff)
        else:
            structChar = TIFFtypeStructCharacter(self.type)

            padding = 4 - self.nValueBytes()
            if padding > 0:
                tiffyLog(TIFFY_LOGLEVEL_DEBUG, "Packing for tag {tag} (type {typ}, count {cnt}) at pos {pos}: {endian}{count}{structchar}{padding}x {vals}".format(tag=self.tag, typ=self.type, cnt=self.count, pos=f.tell(), endian=byteOrder, count=self.count, structchar=structChar, padding=padding, vals=self.values))
                buff = struct.pack("{endian}{count}{structchar}{padding}x".format(endian=byteOrder, count=self.count, structchar=structChar, padding=padding), *self.values)
                f.write(buff)
            else:
                tiffyLog(TIFFY_LOGLEVEL_DEBUG, "Packing for tag {tag} (type {typ}, count {cnt}) at pos {pos}: {endian}{count}{structchar} {vals}".format(tag=self.tag, typ=self.type, cnt=self.count, pos=f.tell(), endian=byteOrder, count=self.count, structchar=structChar, vals=self.values))
                buff = struct.pack("{endian}{count}{structchar}".format(endian=byteOrder, count=self.count, structchar=structChar), *self.values)
                f.write(buff)

    def writeExtraValues(self, f, byteOrder):
        structChar = TIFFtypeStructCharacter(self.type)
        if self.nValueBytes() > 4:
            tiffyLog(TIFFY_LOGLEVEL_DEBUG, "Packing for tag {tag} (type {typ}, count {cnt}) at pos {pos}: {endian}{count}{structchar} {vals}".format(tag=self.tag, typ=self.type, cnt=self.count, pos=f.tell(), endian=byteOrder, count=self.count, structchar=structChar, vals=self.values))

            buff = struct.pack("{endian}{count}{structchar}".format(endian=byteOrder, count=self.count, structchar=structChar), *self.values)
            f.write(buff)


class ifd:
    def __init__(self):
        self.ifdPos = None
        self.fieldBytes = []
        self.fields = []
        self.fieldCount = 0
        self.nextIfdPos = 0

    def addEntry(self, entry):
        self.fields.append(entry)

    def read(self, ifdPos, f, byteOrder):
        self.ifdPos = ifdPos

        tiffyLog(TIFFY_LOGLEVEL_INFO, " -- new IFD.")

        f.seek(ifdPos)
        buff = f.read(2)
        offset = ifdPos + 2
        (self.fieldCount,) = struct.unpack("{endian}H".format(endian=byteOrder), buff)
        
        for i in range(self.fieldCount):
            entry = ifdEntry()
            entry.read(offset, f, byteOrder)
            self.fields.append(entry)
            offset += 12

        f.seek(offset)
        buff = f.read(4)
        (self.nextIfdPos,) = struct.unpack("{endian}L".format(endian=byteOrder), buff)

    def sizeInBytes(self):
        size = 2 + 4  # Number of entries (2 bytes) + offset to next IFD (4 bytes)
        for field in self.fields:
            size += field.sizeInBytes()

        return size

    def prepareDataOffsets(self, offset):
        self.ifdPos = offset

        offset += 2  # number of entries (2 bytes)
        extraDataOffset = offset + 12*len(self.fields) + 4  # offset for each entry (12 bytes) + pointer to next IFD (4 bytes)

        for field in self.fields:
            offset, extraDataOffset = field.prepareDataOffsets(offset, extraDataOffset)

        return offset

    def printOffset(self):
        tiffyLog(TIFFY_LOGLEVEL_INFO, "+ IFD                   {offset} -> {nextIFD}".format(offset=self.ifdPos, nextIFD=self.nextIfdPos))

        for entry in self.fields:
            entry.printOffset()

        for entry in self.fields:
            entry.printExtraDataOffset()

    def write(self, f, byteOrder):
        buff = struct.pack("{endian}H".format(endian=byteOrder), len(self.fields))
        f.write(buff)

        tiffyLog(TIFFY_LOGLEVEL_DEBUG, "Written IFD. Now at position {}".format(f.tell()))

        # IFD entries:
        for entry in self.fields:
            entry.write(f, byteOrder)

        buff = struct.pack("{endian}L".format(endian=byteOrder), self.nextIfdPos)
        f.write(buff)

        # extra data
        for entry in self.fields:
            entry.writeExtraValues(f, byteOrder)


class tiffSubfile:
    def __init__(self):
        self.px = []      # Array of numpy arrays to store image data (i.e. pixel values). One for each component/channel. Shape: (nChannel, nRows, nCols)
        self.imageDataOffset = 0  # Location of beginning of image data
        self.sampleFormat = [TIFF_SAMPLEFORMAT_UINT]  # Standard: uint

        self.filename = None
        self.byteOrder = "<"  # little endian
        self.ifd = None

        self.photometricInterpretation = TIFF_BLACK_IS_ZERO
        self.compression = TIFF_NO_COMPRESSION
        self.predictor = TIFF_NO_PREDICTOR

        self.rows = 0
        self.cols = 0

        self.resolutionUnit = TIFF_RES_INCH
        self.resX = 0
        self.resY = 0

        self.rowsPerStrip = 0
        self.stripOffsets = []
        self.stripByteCounts = []

        self.bitsPerSample = (1,)        # Bilevel images do not define this, so make 1 bit/sample the default.
        self.samplesPerPixel = 1
        self.planarConfig = TIFF_CHUNKY  # Chunky (RGBRGBRGB...)
        self.orientation = 1
        self.colorMap = None             # For palette-color images

    def reset(self):
        self.__init__()

    def set(self, imageData, resX=0, resY=0):
        if len(imageData) > 0:
            self.reset()
            self.px = imageData

            shp = numpy.shape(self.px)  # Shape must be 3-component tuple: (nChannels, height, width)
            if len(shp) == 3:
                self.samplesPerPixel = shp[0]
                self.rows = shp[1]
                self.cols = shp[2]

                self.rowsPerStrip = self.rows
                # self.stripByteCounts and self.stripOffsets will be set later by self.prepareDataOffsets()

                # Byte order
                self.byteOrder = getByteOrder(self.px)

                # Sample format
                if numpy.issubdtype(self.px.dtype, numpy.signedinteger):
                    self.sampleFormat = [TIFF_SAMPLEFORMAT_INT] * self.nChannels()
                elif numpy.issubdtype(self.px.dtype, numpy.unsignedinteger):
                    self.sampleFormat = [TIFF_SAMPLEFORMAT_UINT] * self.nChannels()
                elif numpy.issubdtype(self.px.dtype, numpy.floating):
                    self.sampleFormat = [TIFF_SAMPLEFORMAT_IEEEFP] * self.nChannels()
                else:
                    raise Exception("Unsupported dtype ({}) of provided image data. Must be an integer or floating point type.".format(numpy.dtype(self.px)))

                self.bitsPerSample  = [self.px.dtype.itemsize*8] * self.samplesPerPixel
                tiffyLog(TIFFY_LOGLEVEL_DEBUG, "Setting bits per sample: {}".format(self.bitsPerSample))
                tiffyLog(TIFFY_LOGLEVEL_DEBUG, "Setting sample format:   {}".format(self.sampleFormat))

                self.resolutionUnit = TIFF_RES_NONE
                self.resX = resX
                self.resY = resY

                self.compression    = TIFF_NO_COMPRESSION

                
                self.photometricInterpretation = TIFF_BLACK_IS_ZERO
                if self.nChannels() == 3:
                    self.photometricInterpretation = TIFF_RGB

                #if self.samplesPerPixel > 1:
                #    self.planarConfig   = TIFF_PLANAR
                #else:
                #    self.planarConfig   = TIFF_CHUNKY
                self.planarConfig   = TIFF_CHUNKY

            else:
                raise Exception("Error setting image data. Please provide a numpy array of shape (nChannels, nRows, nColumns).")

    def addIFDentry_shortOrLong(self, tag, values):
        # find the right integer type for the provided values:
        m = max(values)
        typ = TIFF_LONG
        if m < (2**16):
            typ = TIFF_SHORT

        self.addIFDentry(tag, typ, values)

    def addIFDentry(self, tag, typ, values):
        entry = ifdEntry()
        entry.set(tag, typ, values)
        self.ifd.addEntry(entry)

    def setupIFD(self):
        self.ifd = None
        self.ifd = ifd()

        self.addIFDentry_shortOrLong(256, (self.cols, ))                        # n columns
        self.addIFDentry_shortOrLong(257, (self.rows, ))                        # n rows
        self.addIFDentry(258, TIFF_SHORT, self.bitsPerSample)                   # bits per sample, already an array
        self.addIFDentry(259, TIFF_SHORT, (TIFF_NO_COMPRESSION, ))               # compression
        self.addIFDentry(262, TIFF_SHORT, (self.photometricInterpretation, ))   # photometric interpretation
        self.addIFDentry(266, TIFF_SHORT, (TIFF_MSB2LSB, ))                      # fill order 
        self.addIFDentry_shortOrLong(273, (self.imageDataOffset, ))             # offset location of strip
        self.addIFDentry(274, TIFF_SHORT, (self.orientation, ))                 # orientation
        self.addIFDentry(277, TIFF_SHORT, (self.samplesPerPixel, ))             # samples per pixel
        self.addIFDentry_shortOrLong(278, (self.rowsPerStrip, ))                # rows per strip
        self.addIFDentry_shortOrLong(279, (self.dataSizeInBytes(), ))            # strip byte counts
        # insert resolution entries here later...
        self.addIFDentry(284, TIFF_SHORT, (self.planarConfig, ))                # planar configuration
       
        #self.addIFDentry(296, TIFF_SHORT, (self.resolutionUnit, ))             # resolution unit
        
        if self.predictor != TIFF_NO_PREDICTOR:
            self.addIFDentry(317, TIFF_SHORT, (self.predictor, ))                   # predictor
        # color map
        
        self.addIFDentry(339, TIFF_SHORT, self.sampleFormat)                    # sample format, already an array

    def dataSizeInBytes(self):
        """ Size of pixel data (in bytes) """
        size = 0
        for bitsPerSample in self.bitsPerSample:
            s = int(bitsPerSample * self.nPixels())
            size += s

        return size/8

    def sizeInBytes(self):
        size = 0
        if self.ifd != None:
            size += self.ifd.sizeInBytes()

        size += self.dataSizeInBytes()

        return size

    def prepareDataOffsets(self, offset):
        if self.ifd != None:
            self.ifd.prepareDataOffsets(offset)

            offset += self.ifd.sizeInBytes()
            self.imageDataOffset = offset  # image data starts here

            self.stripOffsets    = (self.imageDataOffset, )
            self.stripByteCounts = self.dataSizeInBytes()

            offset += self.dataSizeInBytes()
            return offset
        else:
            return 0

    def printOffset(self):
        if self.ifd != None:
            self.ifd.printOffset()

        tiffyLog(TIFFY_LOGLEVEL_INFO, "  Image Data            {offset}".format(offset=self.imageDataOffset))

    def readMetaInformation(self, filename, byteOrder, imgFileDirectory):
        self.filename = filename
        self.byteOrder = byteOrder
        self.ifd = imgFileDirectory

        # Interpret IFD fields based on their TIFF tags:
        for field in self.ifd.fields:
            if field.tag == 256:  # Number of columns
                self.cols = field.getValue()
            elif field.tag == 257:  # Number of rows
                self.rows = field.getValue()
            elif field.tag == 258:  # Bits per sample
                self.bitsPerSample = field.values
            elif field.tag == 259:  # Compression
                self.compression = field.getValue()
                if self.compression == 0:
                    self.compression = TIFF_NO_COMPRESSION
            elif field.tag == 262:    # Photometric interpretation
                self.photometricInterpretation = field.getValue()
            elif field.tag == 273:  # Strip offsets
                self.stripOffsets = field.values
            elif field.tag == 274:  # Orientation
                self.orientation = field.getValue()
            elif field.tag == 277:  # Samples per pixel
                self.samplesPerPixel = field.getValue()
            elif field.tag == 278:  # Rows per strip
                self.rowsPerStrip = field.getValue()
            elif field.tag == 279:  # Strip byte counts
                self.stripByteCounts = field.values
            elif field.tag == 282:  # x resolution
                self.resX = field.getValue()
            elif field.tag == 283:  # y resolution
                self.resY = field.getValue()
            elif field.tag == 284:  # Planar configuration (order of pixel components)
                self.planarConfig = field.getValue()
            elif field.tag == 296:  # Resolution unit
                self.resolutionUnit = field.getValue()
            elif field.tag == 317:  # Predictor
                self.predictor = field.getValue()
            elif field.tag == 320:  # Color map
                pass # implement later...
            elif field.tag == 339:  # Sample format
                self.sampleFormat = field.values            

    def nChannels(self):
        return self.samplesPerPixel

    def nPixels(self):
        return self.rows*self.cols

    def bitsPerPixel(self):
        bits = 0
        for b in self.bitsPerSample:
            bits += b

        return bits

    def isPlanar(self):
        return (self.samplesPerPixel == 1 or self.planarConfig == TIFF_PLANAR)

    def isSet(self):
        """ Check if image has a valid width and height. """
        if(self.getHeight() > 0):
            if(self.getWidth() > 0):
                return True

        return False

    def getWidth(self):
        return self.cols

    def getHeight(self):
        return self.rows

    def rot90(self):
        if self.isSet():
            self.px = numpy.rot90(self.px, k=1, axes=(1,2))
            self.cols, self.rows = self.rows, self.cols

    def rot180(self):
        if self.isSet():
            self.px = numpy.rot90(self.px, k=2, axes=(1,2))

    def rot270(self):
        if self.isSet():
            self.px = numpy.rot90(self.px, k=-1, axes=(1,2))
            self.cols, self.rows = self.rows, self.cols

    def rotate(self, rotation):
        if rotation == "90":
            self.rot90()
        elif rotation == "180":
            self.rot180()
        elif rotation == "270":
            self.rot270()

    def flipHorizontal(self):
        if self.isSet():
            for i in range(self.nChannels()):
                self.px[i] = numpy.fliplr(self.px[i])

    def flipVertical(self):
        if self.isSet():
            for i in range(self.nChannels()):
                self.px[i] = numpy.flipud(self.px[i])

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

    def numpyDatatype(self, sampleFormat, bitsPerSample):
        formatString = self.byteOrder
        if sampleFormat == TIFF_SAMPLEFORMAT_UINT:  # unsigned
            if bitsPerSample == 8:
                formatString = "u1"   # unsigned char (1 byte)
            elif bitsPerSample == 16:
                formatString += "u2"   # unsigned short (2 bytes)
            elif bitsPerSample == 32:
                formatString += "u4"   # unsigned long (4 bytes)
            elif bitsPerSample == 64:
                formatString += "u8"   # unsigned long long (8 bytes)
        elif sampleFormat == TIFF_SAMPLEFORMAT_INT:  # signed
            if bitsPerSample == 8:
                formatString = "i1"   # signed char (1 byte)
            elif bitsPerSample == 16:
                formatString += "i2"   # signed short (2 bytes)
            elif bitsPerSample == 32:
                formatString += "i4"   # signed long (4 bytes)
            elif bitsPerSample == 64:
                formatString += "i8"   # signed long long (8 bytes)
        elif sampleFormat == TIFF_SAMPLEFORMAT_IEEEFP:
            if bitsPerSample == 16:
                formatString += "f2"   # float 16 bit
            elif bitsPerSample == 32:
                formatString += "f4"   # float 32 bit
            elif bitsPerSample == 64:
                formatString += "f8"   # double 64 bit

        if len(formatString) >= 1:
            return numpy.dtype(formatString)

        raise Exception("Unsupported data type: {bps} bits per sample for TIFF sample format {sf}.".format(bps=bitsPerSample, sf=sampleFormat))

    def structDataTypeString(self, sampleFormat, bitsPerSample):
        formatString = self.byteOrder
        if sampleFormat == TIFF_SAMPLEFORMAT_UINT:  # unsigned
            if bitsPerSample == 8:
                formatString = "B"   # unsigned char (1 byte)
            elif bitsPerSample == 16:
                formatString += "H"   # unsigned short (2 bytes)
            elif bitsPerSample == 32:
                formatString += "L"   # unsigned long (4 bytes)
            elif bitsPerSample == 64:
                formatString += "Q"   # unsigned long long (8 bytes)
        elif sampleFormat == TIFF_SAMPLEFORMAT_INT:  # signed
            if bitsPerSample == 8:
                formatString = "b"   # signed char (1 byte)
            elif bitsPerSample == 16:
                formatString += "h"   # signed short (2 bytes)
            elif bitsPerSample == 32:
                formatString += "l"   # signed long (4 bytes)
            elif bitsPerSample == 64:
                formatString += "q"   # signed long long (8 bytes)
        elif sampleFormat == TIFF_SAMPLEFORMAT_IEEEFP:
            if bitsPerSample == 16:
                formatString += "e"   # float 16 bit
            elif bitsPerSample == 32:
                formatString += "f"   # float 32 bit
            elif bitsPerSample == 64:
                formatString += "d"   # double 64 bit

        if len(formatString) >= 1:
            return formatString

        raise Exception("Unsupported data type: {bps} bits per sample for TIFF sample format {sf}.".format(bps=bitsPerSample, sf=sampleFormat))

    def importFromUncompressedBuffer(self, pixelOffset, nPixels, datatype, channel, buff):
        if self.isPlanar():  # Buffer contains data just for one channel.
            self.px[channel][int(pixelOffset):int(pixelOffset+nPixels)] = numpy.frombuffer(buff, dtype=datatype)
        else: # CHUNKY configuration. channel is irrelevant here (and wrong.)
            bitsPerPixel = self.bitsPerPixel()
            bytesPerPixel = int(bitsPerPixel / 8)

            bytesPerSample = [x/8 for x in self.bitsPerSample]

            buff1d = numpy.frombuffer(buff, dtype=datatype)
            self.px[int(pixelOffset):int(pixelOffset+nPixels)] = numpy.reshape(buff1d, (nPixels, self.nChannels()))

            """
            for pixel in range(nPixels):
                channelOffset = 0

                

                for c in range(self.samplesPerPixel):
                    # Current index in bytes
                    idxStart = int((pixel + pixelOffset)*bytesPerPixel + channelOffset/8)
                    idxStop  = int(idxStart + bytesPerSample[c])

                    #print("Start: {}, Stop: {}".format(idxStart, idxStop))

                    self.px[c][pixel] = numpy.frombuffer(buff[idxStart:idxStop], dtype=dt)

                    channelOffset += self.bitsPerSample[c]
            """

    def imageData(self, obeyOrientation=True):
        # Read data into byte buffer:
        if len(self.stripOffsets) > 0:
            if len(self.stripOffsets) == len(self.stripByteCounts):
                if os.path.isfile(self.filename):
                    with open(self.filename, "rb") as f:
                        # Initialize nChannels x nPixels array for each channel:

                        # Check if bits per sample is the same for each channel:
                        if self.bitsPerSample.count(self.bitsPerSample[0]) == len(self.bitsPerSample) and self.sampleFormat.count(self.sampleFormat[0]) == len(self.sampleFormat):
                            sampleFormat  = self.sampleFormat[0]
                            bitsPerSample = self.bitsPerSample[0]

                            dt = self.numpyDatatype(sampleFormat, bitsPerSample)

                            if self.isPlanar():
                                self.px = numpy.zeros((self.nChannels(), self.nPixels()), dtype=dt)
                            else:
                                # Make array the shape of a chunky tiff configuration and reshape later...
                                self.px = numpy.zeros((self.nPixels(), self.nChannels()), dtype=dt)

                            nStrips = len(self.stripOffsets)

                            if self.compression == TIFF_NO_COMPRESSION or self.compression == TIFF_LZW_COMPRESSION:
                                pixel = 0
                                bitsPerPixel = self.bitsPerPixel()

                                c = 0  # channel id. Only necessary for PLANAR configuration.
                                nStripsPerChannel = int(nStrips / self.samplesPerPixel)

                                for i in range(nStrips):
                                    if self.samplesPerPixel > 1:
                                        if self.planarConfig == TIFF_PLANAR:
                                            # Import next channel once all strips for one component channel have been imported.
                                            if (i % nStripsPerChannel) == 0:
                                                c += 1

                                    offset = self.stripOffsets[i]
                                    byteCount = self.stripByteCounts[i]  # byte count after compression

                                    #print("Strip #{}/{}, Byte Offset: {}, Byte Count: {}".format(i, nStrips, offset, byteCount))

                                    f.seek(offset)
                                    buff = f.read(byteCount)
                                    if self.compression == TIFF_LZW_COMPRESSION:
                                        compressed = lzwData()
                                        compressed.setCompressed(buff)

                                        #print("Decompressing LZW for strip {i}/{n}...".format(i=i, n=nStrips))
                                        compressed.decompress()
                                        buff = bytes(compressed.uncompressed)

                                    byteCount = len(buff)
                                    bitCount = 8*byteCount
                                    if not self.isPlanar(): # Standard is chunky, also for only 1 sample/pixel.
                                        pixelsInStrip = int(bitCount / bitsPerPixel)
                                    else:
                                        pixelsInStrip = int(bitCount / self.bitsPerSample[c])

                                    self.importFromUncompressedBuffer(
                                        pixelOffset = pixel,
                                        nPixels = pixelsInStrip,
                                        datatype = dt,
                                        channel = c,
                                        buff = buff)

                                    pixel += pixelsInStrip

                                tiffyLog(TIFFY_LOGLEVEL_INFO, "{} strips imported.".format(nStrips))
                            else:
                                raise Exception("TIFF: Compression scheme {} not supported.".format(self.compression))

                            f.close()

                            if not(self.isPlanar()):
                                # Chunky mode. Swap axes from (pixels, channels) to (channels, pixels):
                                self.px = numpy.swapaxes(self.px, 0, 1)

                            # Reshape components into 2D arrays:
                            if len(self.px) > 0:
                                self.px = numpy.reshape(self.px, (self.nChannels(), self.rows, self.cols))

                            # Convert horizontal differences to absolute values:
                            tiffyLog(TIFFY_LOGLEVEL_INFO, "Applying horizontal differencing...")
                            if self.predictor == TIFF_HORIZONTAL_DIFFERENCING:
                                for col in range(1, self.cols):
                                    self.px[...,col] = self.px[...,col-1] + self.px[...,col]

                            tiffyLog(TIFFY_LOGLEVEL_INFO, "Orientation is: {}".format(self.orientation))
                            if obeyOrientation:
                                # Rotate back to orientation 1 ((0,0) is upper left)
                                if self.orientation == 2:     # (col0, row0) is (top, right)
                                    self.flip(horizontal=True, vertical=False)
                                elif self.orientation == 3:   # (col0, row0) is (bottom, right)
                                    self.flip(horizontal=True, vertical=True)
                                elif self.orientation == 4:   # (col0, row0) is bottom left
                                    self.flip(horizontal=False, vertical=True)
                                elif self.orientation == 5:   # (col0, row0) is (left, top)
                                    self.rotate("90")
                                    self.flip(horizontal=False, vertical=True)
                                elif self.orientation == 6:   # (col0, row0) is (right, top)
                                    self.rotate("270")
                                elif self.orientation == 7:   # (col0, row0) is (right, bottom)
                                    self.rotate("90")
                                    self.flip(horizontal=True, vertical=False)
                                elif self.orientation == 8:   # (col0, row0) is (left, bottom)
                                    self.rotate("90")

                                self.orientation = 1

                            return self.px
                        else:
                            f.close()
                            raise Exception("Unsupported TIFF format: all channels must have the same type and size. Sample format: {}, Bits per sample: {}".format(self.sampleFormat, self.bitsPerSample))                       

                raise Exception("File not available: {}".format(self.filename))
            else:
                raise Exception("Number of strip offsets ({nOffsets}) does not match number of strip byte counts ({nByteCounts}).".format(nOffsets=len(self.stripOffsets), nByteCounts=len(self.stripByteCounts)))
        else:
            raise Exception("No data strips found for requested subfile in {filename}.".format(filename=self.filename))

    def write(self, f, byteOrder):
        """ Expects an open, writable file pointer f. """
        self.ifd.write(f, byteOrder)

        dataByteOrder = getByteOrder(self.px)
        if dataByteOrder != byteOrder:
            self.px.byteswap(inplace=True)

        # Write image data:
        if self.samplesPerPixel == 1:
            f.write(self.px)
        else:
            chunkyBytes = numpy.swapaxes(self.px, 0, 2)  # channel <-> cols  --> (cols, rows, channel)
            chunkyBytes = numpy.swapaxes(chunkyBytes, 0, 1)  # cols <-> rows  --> (rows, cols, channel)
            chunkyBytes.tofile(f, "")

            """
            dataString = []
            for c in range(self.nChannels()):
                dataString.append(self.structDataTypeString(self.sampleFormat[c], self.bitsPerSample[c]))

            for y in range(self.rows):
                for x in range(self.cols):
                    # Chunky style.
                    for c in range(self.nChannels()):
                        buff = struct.pack("{endian}{ds}".format(endian=byteOrder, ds=dataString[c]), self.px[c][y][x])
                        f.write(buff)
            """

class tiff:
    def __init__(self):
        self.filename = None
        self.byteOrder = "<"
        self.subfiles = []

    def reset(self):
        self.__init__()

    def read(self, filename=None):
        self.filename = filename

        if self.filename != None:
            if os.path.isfile(self.filename):
                filesize = os.path.getsize(self.filename)
                tiffyLog(TIFFY_LOGLEVEL_DEBUG, "Size of {filename}: {size} Bytes.".format(filename=self.filename, size=filesize))

                if filesize > 2:
                    with open(self.filename, "rb") as f:
                        buff = f.read(2)
                        (tiffformat,) = struct.unpack("<H", buff)  # Just read as little endian. It's symmetric, it doesn't matter.

                        self.byteOrder = "<"   # little endian
                        if tiffformat == 0x4949:
                            tiffyLog(TIFFY_LOGLEVEL_DEBUG, "II format. Little endian.")
                            self.byteOrder = "<"
                        elif tiffformat == 0x4d4d:
                            tiffyLog(TIFFY_LOGLEVEL_DEBUG, "MM format. Big endian.")
                            self.byteOrder = ">"
                        else:
                            f.close()
                            raise Exception("Invalid byte order for first two bytes in TIFF header. Must be 0x4949 (II) or 0x4d4d (MM).")

                        buff = f.read(2)
                        (magicByte,) = struct.unpack("{endian}H".format(endian=self.byteOrder), buff)
                        if magicByte == 42:
                            tiffyLog(TIFFY_LOGLEVEL_DEBUG, "Magic Byte: {}".format(magicByte))
                        else:
                            raise Exception("TIFF magic byte is not 42: {}".format(self.filename))


                        # Get location of first image file directory (IFD):
                        buff = f.read(4)
                        (ifdPos,) = struct.unpack("{endian}L".format(endian=self.byteOrder), buff)

                        ifd0 = ifd()
                        ifd0.read(ifdPos, f, self.byteOrder)

                        self.subfiles = []
                        subfile0 = tiffSubfile()
                        subfile0.readMetaInformation(self.filename, self.byteOrder, ifd0)
                        self.subfiles.append(subfile0)

                        nextIfdPos = ifd0.nextIfdPos
                        while nextIfdPos != 0:
                            ifdn = ifd()
                            ifdn.read(nextIfdPos, f, self.byteOrder)
                            subfilen = tiffSubfile()
                            subfilen.readMetaInformation(self.filename, self.byteOrder, ifdn)
                            self.subfiles.append(subfilen)
                            nextIfdPos = ifdn.nextIfdPos

                        f.close()
                else:
                    raise Exception("TIFF file does not contain any data: {}".format(self.filename))
            else:
                raise Exception("File not found: {}".format(self.filename))

    def imageData(self, subfile=0, channel=None, obeyOrientation=True):
        if subfile < len(self.subfiles):
            data = self.subfiles[subfile].imageData()
            if channel is None:
                return data
            else:
                return data[0]
        else:
            raise Exception("Subfile id not available: {}".format(subfile))

    def getNrSubfiles(self):
        """ Return number of subfiles. """
        return len(self.subfiles)

    def getWidth(self, subfile=0):
        if subfile < len(self.subfiles):
            return self.subfiles[subfile].getWidth()

    def getHeight(self, subfile=0):
        if subfile < len(self.subfiles):
            return self.subfiles[subfile].getHeight()

    def sizeInBytes(self):
        size = 8  # TIFF header

        # Add size for each sub file:
        for sub in self.subfiles:
            size += sub.sizeInBytes()

        return size

    def prepareDataOffsets(self):
        offset = 8  # TIFF header
        for sub in self.subfiles:
            offset = sub.prepareDataOffsets(offset)

    def printOffset(self):
        self.prepareDataOffsets()
        tiffyLog(TIFFY_LOGLEVEL_INFO, "TIFF HEADER             0")
        for sub in self.subfiles:
            sub.printOffset()

    def set(self, imageData, resX=0, resY=0):
        self.reset()
        self.addImgData(imageData, resX, resY)

    def addImgData(self, imageData, resX=0, resY=0):
        if len(numpy.shape(imageData)) == 2:
            # We have to add another dimension for the channel...
            imageData = numpy.expand_dims(a=imageData, axis=0)

        if len(numpy.shape(imageData)) == 3:
            img = tiffSubfile()
            img.set(imageData, resX, resY)
            self.subfiles.append(img)

            # IFDs must be setup twice to ensure correct data offset pointers.
            for sub in self.subfiles:
                sub.setupIFD()
            self.prepareDataOffsets()
            for sub in self.subfiles:
                sub.setupIFD()
        else:
            raise Exception("TIFF: adding image data failed. Image data must be a numpy array of shape (height, width) or (channels, height, width).")

    def save(self, filename, endian="little"):
        byteOrder = "<" # little endian is default
        if endian == "big":
            byteOrder = ">"

        with open(filename, 'wb') as f:
            # Write the TIFF header.

            # Endian:
            if byteOrder == ">":
                buff = struct.pack("{endian}H".format(endian=byteOrder), 0x4d4d)
            else:
                buff = struct.pack("{endian}H".format(endian=byteOrder), 0x4949)
            f.write(buff)

            # Magic number:
            buff = struct.pack("{endian}H".format(endian=byteOrder), 42)
            f.write(buff)

            # The IFD0 always follows the header in our case:
            buff = struct.pack("{endian}L".format(endian=byteOrder), 8)
            f.write(buff)

            for sub in self.subfiles:
                sub.write(f, byteOrder)
           
            f.close()
