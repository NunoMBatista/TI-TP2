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

	def readDinamicBlock (self):
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

	def storeLITLENcodeLens(self, HLIT, CLENTree):
		'''Takes the code lengths huffmantree and stores the 
		HLIT + 257 lit/comp code lengths accordingly'''

		# Array where the code lengths will be stored 
		LITLENcodeLens = [] 
		#inarray = 0
		prevCode = -1
		while (len(LITLENcodeLens) < HLIT + 257):
			CLENTree.resetCurNode()
			found = False
			while(not found):
				curBit = self.readBits(1)
				code = CLENTree.nextNode(str(curBit))
				if(code != -1 and code != -2):
					found = True
	   
			if(code == 18):
				ammount = self.readBits(7)
				LITLENcodeLens += [0]*(11 + ammount)
			if(code == 17):
				ammount = self.readBits(3)
				LITLENcodeLens += [0]*(3 + ammount)
			if(code == 16):
				ammount = self.readBits(2)
				LITLENcodeLens += [prevCode]*(3 + ammount)
			elif(code >= 0 and code <= 15):
				LITLENcodeLens += [code]
				prevCode = code
    
		return LITLENcodeLens

	def storeDISTcodeLens(self, HDIST, CLENTree):
		'''Takes the code lengths huffmantree and stores the 
		HDIST distance code lengths accordingly'''

		# Array where the code lengths will be stored 
		DISTcodeLens = [] 
		#inarray = 0
		prevCode = -1
		while (len(DISTcodeLens) < HDIST + 1):
			CLENTree.resetCurNode()
			found = False
			while(not found):
				curBit = self.readBits(1)
				code = CLENTree.nextNode(str(curBit))
				if(code != -1 and code != -2):
					found = True
     
			if(code == 18):
				ammount = self.readBits(7)
				DISTcodeLens += [0]*(11 + ammount)
			if(code == 17):
				ammount = self.readBits(3)
				DISTcodeLens += [0]*(3 + ammount)
			if(code == 16):
				ammount = self.readBits(2)
				DISTcodeLens += [prevCode]*(3 + ammount)
			elif(code >= 0 and code <= 15):
				DISTcodeLens += [code]
				prevCode = code

		return DISTcodeLens

	def decompressLZ77(self, HuffmanTreeLITLEN, HuffmanTreeDIST):
		ExtraLITLENBits = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5]
		ExtraLITLENLens = [11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227]
		ExtraDISTBits = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13]		
		ExtraDISTLens = [5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577]

		codeLITLEN = -1
		output = []
		while(codeLITLEN != 256):
			HuffmanTreeLITLEN.resetCurNode()
			foundLITLEN = False
			distFound = True

			while(not foundLITLEN):
				curBit = str(self.readBits(1))
				codeLITLEN = HuffmanTreeLITLEN.nextNode(curBit)
    
				if (codeLITLEN != -1 and codeLITLEN != -2):
					foundLITLEN = True
					if(codeLITLEN < 256):
						output += [codeLITLEN]
					if(codeLITLEN > 256):
						distFound = False
						if(codeLITLEN < 265):
							length = codeLITLEN - 257 + 3
						else:
							dif = codeLITLEN - 265
							readExtra = ExtraLITLENBits[dif]
							lenExtra = ExtraLITLENLens[dif]
							length = lenExtra + self.readBits(readExtra)

						HuffmanTreeDIST.resetCurNode()
						while(not distFound):
							distBit = str(self.readBits(1))
							codeDIST = HuffmanTreeDIST.nextNode(distBit)
       
							if(codeDIST != -1 and codeDIST != -2):
								distFound = True
								if(codeDIST < 4):
									distance = codeDIST + 1

								else:
									dif = codeDIST - 4
									readExtra = ExtraDISTBits[dif]
									distExtra = ExtraDISTLens[dif]
									distance = distExtra + self.readBits(readExtra)
								
								#outputCopy = output[-distance+2:-distance+2+length]

								for i in range(length):
									output.append(output[-distance])
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
		while not BFINAL == 1:	
      
			BFINAL = self.readBits(1)
			
			BTYPE = self.readBits(2)					
			if BTYPE != 2:
				print('Error: Block %d not coded with Huffman Dynamic coding' % (numBlocks+1))
				return
			
			# if BTYPE == 10 in base 2 -> read the dinamic Huffman compression format 
			if BTYPE == int('10', 2):		
				# HLIT: # of literal/length  codes
				# HDIST: # of distance codes 
				# HCLEN: # of code length codes
				HLIT, HDIST, HCLEN = self.readDinamicBlock()
				
				# print("HLIT:", bin(HLIT)[2:], "=", HLIT, 
          		# 	"\nHDIST:", bin(HDIST)[2:], "=", HDIST,
             	# 	"\nHCLEN:", bin(HCLEN)[2:], "=", HCLEN)

				CLENcodeLens = self.storeCLENLengths(HCLEN)   
				#print("Code Lengths of indices i from the code length tree:", CLENcodeLens)

				HuffmanTreeCLENs = self.createHuffmanFromLens(CLENcodeLens, verbose=False)
				LITLENcodeLens = self.storeLITLENcodeLens(HLIT, HuffmanTreeCLENs)
				
				HuffmanTreeLITLEN = self.createHuffmanFromLens(LITLENcodeLens, verbose=False)

				# for i in range(0, len(LITLENcodeLens)):
				# 	print(i, ":", LITLENcodeLens[i])
     
				DISTcodeLens = self.storeDISTcodeLens(HDIST, HuffmanTreeCLENs)
				# for i in range(0, len(DISTcodeLens)):
				# 	print(i, ":", DISTcodeLens[i])
     
				HuffmanTreeDIST = self.createHuffmanFromLens(DISTcodeLens, verbose=False)
    
				output = self.decompressLZ77(HuffmanTreeLITLEN, HuffmanTreeDIST)

				# for i in output:
				# 	print(i, end="")
			# update number of blocks read
			numBlocks += 1
   
		f = open(self.gzh.fName, 'wb')
		f.write(bytes(output))
		f.close()

		# close file			
		
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
	