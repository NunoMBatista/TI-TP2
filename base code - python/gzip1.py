# Author: Marco Simoes
# Adapted from Java's implementation of Rui Pedro Paiva
# Teoria da Informacao, LEI, 2022

import sys
from huffmantree import HuffmanTree


class GZIPHeader:
	''' class for reading and storing GZIP header fields '''

	ID1 = ID2 = CM = FLG = XFL = OS = 0
	MTIME = []
	lenMTIME = 4
	mTime = 0

	# bits 0, 1, 2, 3 and 4, respectively (remaining 3 bits: reserved)
	FLG_FTEXT = FLG_FHCRC = FLG_FEXTRA = FLG_FNAME = FLG_FCOMMENT = 0   
	
	# FLG_FTEXT --> ignored (usually 0)
	# if FLG_FEXTRA == 1
	XLEN, extraField = [], []
	lenXLEN = 2
	
	# if FLG_FNAME == 1
	fName = ''  # ends when a byte with value 0 is read
	
	# if FLG_FCOMMENT == 1
	fComment = ''   # ends when a byte with value 0 is read
		
	# if FLG_HCRC == 1
	HCRC = []
		
		
	
	def read(self, f):
		''' reads and processes the Huffman header from file. Returns 0 if no error, -1 otherwise '''

		# ID 1 and 2: fixed values
		self.ID1 = f.read(1)[0]  
		if self.ID1 != 0x1f: return -1 # error in the header
			
		self.ID2 = f.read(1)[0]
		if self.ID2 != 0x8b: return -1 # error in the header
		
		# CM - Compression Method: must be the value 8 for deflate
		self.CM = f.read(1)[0]
		if self.CM != 0x08: return -1 # error in the header
					
		# Flags
		self.FLG = f.read(1)[0]
		
		# MTIME
		self.MTIME = [0]*self.lenMTIME
		self.mTime = 0
		for i in range(self.lenMTIME):
			self.MTIME[i] = f.read(1)[0]
			self.mTime += self.MTIME[i] << (8 * i) 				
						
		# XFL (not processed...)
		self.XFL = f.read(1)[0]
		
		# OS (not processed...)
		self.OS = f.read(1)[0]
  
		# --- Check Flags
		self.FLG_FTEXT = self.FLG & 0x01
		self.FLG_FHCRC = (self.FLG & 0x02) >> 1
		self.FLG_FEXTRA = (self.FLG & 0x04) >> 2
		self.FLG_FNAME = (self.FLG & 0x08) >> 3
		self.FLG_FCOMMENT = (self.FLG & 0x10) >> 4
					
		# FLG_EXTRA
		if self.FLG_FEXTRA == 1:
			# read 2 bytes XLEN + XLEN bytes de extra field
			# 1st byte: LSB, 2nd: MSB
			self.XLEN = [0]*self.lenXLEN
			self.XLEN[0] = f.read(1)[0]
			self.XLEN[1] = f.read(1)[0]
			self.xlen = self.XLEN[1] << 8 + self.XLEN[0]
			
			# read extraField and ignore its values
			self.extraField = f.read(self.xlen)
		
		def read_str_until_0(f):
			s = ''
			while True:
				c = f.read(1)[0]
				if c == 0: 
					return s
				s += chr(c)
		
		# FLG_FNAME
		if self.FLG_FNAME == 1:
			self.fName = read_str_until_0(f)
		
		# FLG_FCOMMENT
		if self.FLG_FCOMMENT == 1:
			self.fComment = read_str_until_0(f)
		
		# FLG_FHCRC (not processed...)
		if self.FLG_FHCRC == 1:
			self.HCRC = f.read(2)
			
		return 0
			
   
class GZIP:
	''' class for GZIP decompressing file (if compressed with deflate) '''

	gzh = None
	gzFile = ''
	fileSize = origFileSize = -1
	numBlocks = 0
	f = None
	

	bits_buffer = 0
	available_bits = 0		


	def __init__(self, filename):
		self.gzFile = filename
		self.f = open(filename, 'rb')
		self.f.seek(0,2)
		self.fileSize = self.f.tell()
		self.f.seek(0)

	def readDynamicBlock (self):
		'''Interprets Dinamic Huffman compressed blocks'''
  
		HLIT = self.readBits(5)
		HDIST = self.readBits(5)
		HCLEN = self.readBits(4)
  
		return HLIT, HDIST, HCLEN

	def storeCLENLengths(self, HCLEN):
		'''Stores the code lengths for the code lengths alphabet in an array'''
     
		# Order of lengths in which the bits are read
		idxCLENcodeLens = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15]
		CLENcodeLens = [0 for i in range(19)]

		# CLENcodeLens[idx] = N translates to: "the code for idx in the code lengths alphabet has a length of N"
		# if N == 0, that indexes' code length is not used
		for i in range(0, HCLEN+4):
			temp = self.readBits(3)
			CLENcodeLens[idxCLENcodeLens[i]] = temp
		return CLENcodeLens

	def createHuffmanFromLens(self, lenArray, verbose=False):
		'''Takes an array with symbols' Huffman codes' lengths and returns
		a formated Huffman tree with said codes
  
		If verbose==True, it prints the codes as they're added to the tree'''
  
		htr = HuffmanTree()
		# max_len is the code with the largest length 
		max_len = max(lenArray)
		# max_symbol é o maior símbolo a codificar
		max_symbol = len(lenArray)
		
		bl_count = [0 for i in range(max_len+1)]
		# Get array with number of codes with length N (bl_count)
		for N in range(1, max_len+1):
			bl_count[N] += lenArray.count(N)

		# Get first code of each code length 
		code = 0
		next_code = [0 for i in range(max_len+1)]
		for bits in range(1, max_len+1):
			code = (code + bl_count[bits-1]) << 1
			next_code[bits] = code
  
		# Define codes for each symbol in lexicographical order
		for n in range(max_symbol):
			# Length associated with symbol n 
			length = lenArray[n]
			if(length != 0):
				code = bin(next_code[length])[2:]
				# In case there are 0s at the start of the code, we have to add them manualy
				# length-len(code) 0s have to be added
				extension = "0"*(length-len(code)) 
				htr.addNode(extension + code, n, verbose)
				next_code[length] += 1
		
		return htr;

	def storeTreeCodeLens(self, size, CLENTree):
		'''Takes the code lengths huffmantree and stores the code lengths accordingly'''

		# Array where the code lengths will be stored 
		treeCodeLens = [] 
  
		while (len(treeCodeLens) < size):
			# Sets the current node to the root of the tree
			CLENTree.resetCurNode()
			found = False
   
			# While reading, if a leaf hasn't been found, keep searching bit by bit
			while(not found):
				curBit = self.readBits(1)
				code = CLENTree.nextNode(str(curBit))
				if(code != -1 and code != -2):
					found = True

			# SPECIAL CHARACTERS
			# 18 - Reads 7 extra bits 
			# 17 - Reads 3 extra bits
			# 16 - Reads 2 extra bits
			if(code == 18):
				amount = self.readBits(7)
				# According to the 7 bits just read, set the following 11-138 values on the length array to 0 
				treeCodeLens += [0]*(11 + amount)
			if(code == 17):
				amount = self.readBits(3)
				# According to the 3 bits just read, set the following 3-10 values on the length array to 0 
				treeCodeLens += [0]*(3 + amount)
			if(code == 16):
				amount = self.readBits(2)
				# According to the 2 bits just read, set the following 3-6 values on the length array to the latest length read
				treeCodeLens += [prevCode]*(3 + amount)
			elif(code >= 0 and code <= 15):
				# If a special character isn't found, just set the next code length to the value found
				treeCodeLens += [code]
				prevCode = code

		return treeCodeLens

	def decompressLZ77(self, HuffmanTreeLITLEN, HuffmanTreeDIST, output):
     
		# How many bits required to read if length code read is larger than 265
		ExtraLITLENBits = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0]
		# Length required to add if the length code read if larger than 265
		ExtraLITLENLens = [11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258]
	
		# How many bits required to read if the distance code read is larger than 4
		ExtraDISTBits = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13]		
		# Distance required to add if the special character read if larger than 4
		ExtraDISTLens = [5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577]

		codeLITLEN = -1

		# Read from the input stream until 256 is found
		while(codeLITLEN != 256):
			# Resets the current node to the tree's root
			HuffmanTreeLITLEN.resetCurNode()

			foundLITLEN = False
			# A distance is considered "found" by default in case the character read is just a literal
			distFound = True

			# While a literal or length isn't found in the LITLEN tree, keep searching bit by bit
			while(not foundLITLEN):
				curBit = str(self.readBits(1))
				# Updates the current node according the the bit just read
				codeLITLEN = HuffmanTreeLITLEN.nextNode(curBit)
    
				# If a leaf is reached in the LITLEN tree, follow the instructions according to the value found
				if (codeLITLEN != -1 and codeLITLEN != -2):
					foundLITLEN = True
     
					# If the code reached is in the interval [0, 256[, just append the value read corresponding a the literal to the output array
					if(codeLITLEN < 256):
						output += [codeLITLEN]

					# If the code is in the interval [257, 285], it is refering to the length of the string to copy
					if(codeLITLEN > 256):
         
						distFound = False
      
						# if the code is in the interval [257, 265[, sets the length of the string to copy to the code read - 257 + 3
						if(codeLITLEN < 265):
							length = codeLITLEN - 257 + 3
       
						# the codes in the interval [265, 285] are special and require more bits to be read
						else:
							# dif defines the indices in the "Extra array's" to be used 
							dif = codeLITLEN - 265
							# How many extra bits will need to be read
							readExtra = ExtraLITLENBits[dif]
							# How much extra length to add
							lenExtra = ExtraLITLENLens[dif]
							length = lenExtra + self.readBits(readExtra)

						# Resets the curent node in the distance tree to it's root
						HuffmanTreeDIST.resetCurNode()
						# While a distance isn's found in the DIST tree, keep searching bit by bit
						while(not distFound):
							distBit = str(self.readBits(1))
							# Updates the current node according to the bit just read
							codeDIST = HuffmanTreeDIST.nextNode(distBit)

							# If a leaf is reached in the LITLEN tree, follow the instructions according to the value found
							if(codeDIST != -1 and codeDIST != -2):
								distFound = True

								# If the code read is in the interval [0, 4[ define the distance to go back to the code read + 1
								if(codeDIST < 4):
									distance = codeDIST + 1

								# The codes in the interval [4, 29] are special and require more bits to be read
								else:
									# dif defines the indices in the "Extra arrays" to be used
									dif = codeDIST - 4
									readExtra = ExtraDISTBits[dif]
									# How many extra bits need to be read
									distExtra = ExtraDISTLens[dif]
									# How much extra distance to add
									distance = distExtra + self.readBits(readExtra)
								
								# For each one of the range(length) iterations, copy the character at index len(output)-distance to the end of the output array
								for i in range(length):
									output.append(output[-distance])
		#print(output)
		return output
 
	def decompress(self):
		''' main function for decompressing the gzip file with deflate algorithm '''
		numBlocks = 0

		# get original file size: size of file before compression
		origFileSize = self.getOrigFileSize()
		print(origFileSize)
		
		# read GZIP header
		error = self.getHeader()
		if error != 0:
			print('Formato invalido!')
			return
		
		# show filename read from GZIP header
		print(self.gzh.fName)
		
		
		# MAIN LOOP - decode block by block
		BFINAL = 0	

		# Opens the output file in "write" binary mode
		f = open(self.gzh.fName, 'wb')		
  
		output = []
		while not BFINAL == 1:	
			BFINAL = self.readBits(1)
			
			# if BTYPE == 10 in base 2 -> read the dinamic Huffman compression format 
			BTYPE = self.readBits(2)					
			if BTYPE != 2:
				print('Error: Block %d not coded with Huffman Dynamic coding' % (numBlocks+1))
				return
			
			# HLIT: # of literal/length  codes
			# HDIST: # of distance codes 
			# HCLEN: # of code length codes
			HLIT, HDIST, HCLEN = self.readDynamicBlock()
			
			# Store the CLEN tree's code lens in a pre-determined order 
			CLENcodeLens = self.storeCLENLengths(HCLEN)   
			#print("Code Lengths of indices i from the code length tree:", CLENcodeLens)
				
			# Based on the CLEN tree's code lens, define an huffman tree for CLEN
			HuffmanTreeCLENs = self.createHuffmanFromLens(CLENcodeLens, verbose=False)

			# Store the literal and length tree code lens based on the CLEN tree codes
			#LITLENcodeLens = self.storeLITLENcodeLens(HLIT, HuffmanTreeCLENs)
			LITLENcodeLens = self.storeTreeCodeLens(HLIT + 257, HuffmanTreeCLENs)			

			# Define the literal and length huffman tree based on the lengths of it's codes
			HuffmanTreeLITLEN = self.createHuffmanFromLens(LITLENcodeLens, verbose=False)
	
			# Store the distance tree code lens based on the CLEN tree codes
			DISTcodeLens = self.storeTreeCodeLens(HDIST + 1, HuffmanTreeCLENs)

			# Define the distance huffman tree based on the lengths of it's codes
			HuffmanTreeDIST = self.createHuffmanFromLens(DISTcodeLens, verbose=False)

			# Based on the trees defined so far, decompress the data according to the Lempel-Ziv77 algorthm 
			output = self.decompressLZ77(HuffmanTreeLITLEN, HuffmanTreeDIST, output)
   

			# Only the last 32768 characters should be kept in memory
			if(len(output) > 32768):
				# Write every charater that exceeds the 327680 range to the file
				f.write(bytes(output[0 : len(output) - 32768]))
				# Keep the rest in the output array
				output = output[len(output) - 32768 :]

			# update number of blocks read
			numBlocks += 1

		# Write the bytes corresponding to the output array elements
		f.write(bytes(output))
  
		# Close the file
		f.close		

		self.f.close()	
		print("End: %d block(s) analyzed." % numBlocks)
	
	
	def getOrigFileSize(self):
		''' reads file size of original file (before compression) - ISIZE '''
		
		# saves current position of file pointer
		fp = self.f.tell()
		
		# jumps to end-4 position
		self.f.seek(self.fileSize-4)
		
		# reads the last 4 bytes (LITTLE ENDIAN)
		sz = 0
		for i in range(4): 
			sz += self.f.read(1)[0] << (8*i)
		
		# restores file pointer to its original position
		self.f.seek(fp)
		
		return sz		
	

	
	def getHeader(self):  
		''' reads GZIP header'''

		self.gzh = GZIPHeader()
		header_error = self.gzh.read(self.f)
		return header_error
		

	def readBits(self, n, keep=False):
		''' reads n bits from bits_buffer. if keep = True, leaves bits in the buffer for future accesses '''

		while n > self.available_bits:
			self.bits_buffer = self.f.read(1)[0] << self.available_bits | self.bits_buffer
			self.available_bits += 8
		
		mask = (2**n)-1
		value = self.bits_buffer & mask

		if not keep:
			self.bits_buffer >>= n
			self.available_bits -= n

		return value

	

if __name__ == '__main__':

	# gets filename from command line if provided
	fileName = "FAQ.txt.gz"
	if len(sys.argv) > 1:
		fileName = sys.argv[1]			

	# decompress file
	gz = GZIP(fileName)
	gz.decompress()
	
