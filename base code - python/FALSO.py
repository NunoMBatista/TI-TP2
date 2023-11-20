# Author: Marco Simoes
# Adapted from Java's implementation of Rui Pedro Paiva
# Teoria da Informacao, LEI, 2022
import gzip
import sys
import numpy as np
from huffmantree import HuffmanTree
from collections import Counter

###############################################################################
# ARRAYS AUXILIARES À ALÍNEA 7 DO TRABALHO

# NESTE ARRAY GUARDAMOS OS BITS EXTRA A LER NO HLIT
BITSHLIT = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5]
# NESTE ARRAY GUARDAMOS OS TAMANHOS BASE NO HLIT
LENGTHSHLIT = [11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163,
               195, 227]


# NESTE ARRAY GUARDAMOS OS BITS EXTRA A LER NO HDIST
BITSHDIST = [0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11,
             11, 12, 12, 13, 13]
# NESTE ARRAY GUARDAMOS OS TAMANHOS BASE NO HDIST
LENGTHSHDIST = [1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513,
                769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577]

###############################################################################


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
    fComment = ''  # ends when a byte with value 0 is read

    # if FLG_HCRC == 1
    HCRC = []

    def read(self, f):
        ''' reads and processes the Huffman header from file. Returns 0 if no error, -1 otherwise '''

        # ID 1 and 2: fixed values
        self.ID1 = f.read(1)[0]
        if self.ID1 != 0x1f:
            return -1  # error in the header

        self.ID2 = f.read(1)[0]
        if self.ID2 != 0x8b:
            return -1  # error in the header

        # CM - Compression Method: must be the value 8 for deflate
        self.CM = f.read(1)[0]
        if self.CM != 0x08:
            return -1  # error in the header

        # Flags
        self.FLG = f.read(1)[0]

        # MTIME
        self.MTIME = [0] * self.lenMTIME
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
            self.XLEN = [0] * self.lenXLEN
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

###############################################################################
# FUNCOES CRIADAS


# lê e armazena num array os x comprimentos dos codigos (ponto4 -> x=HLIT+257 || ponto5 -> x=HDIST+1)
def defBitABit(self, x, CODIGOS_BINARIO):
    LISTA = []
    # le bit a bit ate encontrar um codigo hufmman em CODIGOS_BINARIO
    strTmp = ""
    contador = 0
    while x != contador:
        strTmp = strTmp + str(self.readBits(1))
        if strTmp in CODIGOS_BINARIO:
            pos = CODIGOS_BINARIO.index(strTmp)

            # se o simbolo for 16,17 ou 18 lê os bits necessários e insere-os no array após o tratamento
            if pos == 16:
                repeticao = self.readBits(2)
                lenAnterior = LISTA[len(LISTA) - 1]
                for j in range(repeticao + 3):
                    LISTA.append(lenAnterior)
                    contador = contador + 1

            elif pos == 17:
                repeticao = self.readBits(3)
                for j in range(repeticao + 3):
                    LISTA.append(0)
                    contador = contador + 1

            elif pos == 18:
                repeticao = self.readBits(7)
                for j in range(repeticao + 11):
                    LISTA.append(0)
                    contador = contador + 1
            # se o simbolo for < 16 insere-os diretamente no array
            else:
                LISTA.append(pos)
                contador = contador + 1

            strTmp = ""

    return LISTA
    pass


# procedimento para obter os codigos de huffman no ponto 3(HCLEN), ponto 4(HLIT) e ponto5(HDIST)
def codigosHuffman(COD, x):
    # com o codigo do Doc1 busca os codigos de huffman(decimal)
    bl_count = Counter(COD)
    bl_count[0] = 0

    next_code = []

    code = 0
    for bits in range(1, max(bl_count) + 1):
        code = (code + bl_count[bits - 1]) << 1
        next_code.append(code)

    CODIGOS = np.zeros(x, dtype=int)

    for n in range(0, x):
        length = COD[n]
        if length != 0:
            CODIGOS[n] = next_code[length - 1]
            next_code[length - 1] = next_code[length - 1] + 1

    # transforma os codigos huffman em binario
    CODIGOS_BINARIO = []
    for i in range(x):
        # obtem o comprimento para determinar se i é simbolo ou trash
        comprimento = COD[i]

        # trash
        if comprimento == 0:
            CODIGOS_BINARIO.append('')

        # simbolo
        else:
            read = CODIGOS[i]

            strDec = ''
            mascara = 1
            # transforma em binario
            for i in range(comprimento):
                strDec = str((read >> i) & mascara) + strDec

            CODIGOS_BINARIO.append(strDec)

    return CODIGOS_BINARIO
    pass

###############################################################################


class GZIP:
    ''' class for GZIP decompressing file (if compressed with deflate) '''

    gzh = None
    gzFile = ''
    fileSize = origFileSize = -1
    numBlocks = 0
    f = None

    bits_buffer = 0
    available_bits = 0

    def __init__(self, filenFALSOame):
        self.gzFile = filename
        self.f = open(filename, 'rb')
        self.f.seek(0, 2)
        self.fileSize = self.f.tell()
        self.f.seek(0)

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
                print('Error: Block %d not coded with Huffman Dynamic coding' % (numBlocks + 1))
                return

            # --- STUDENTS --- ADD CODE HERE

            # -----------------------------------------------------------------
            # SEMANA 1, PONTO 1
            # LER O FORMATO DO BLOCO
            # SEGUNDO A ESPECIFICAO RFC
            # -----------------------------------------------------------------

            HLIT = self.readBits(5) + 257
            HDIST = self.readBits(5) + 1
            HCLEN = self.readBits(4) + 4

            # print("HLIT -> " + str(HLIT))
            # print("HDIST -> " + str(HDIST))
            # print("HCLEN -> " + str(HCLEN))
            # print("TAMANHO HCLEN -> " + str(HCLEN))





            # -----------------------------------------------------------------
            # SEMANA 2, PONTO 2
            # ARMAZENAR OS COMPRIMENTOS DE CODIGO DE FORMA ORDENADA
            # -----------------------------------------------------------------

            HCLEN_ORDENADO = np.zeros(19, dtype=int)

            ordens = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15]
            for i in range(HCLEN):
                HCLEN_ORDENADO[ordens[i]] = self.readBits(3)

            # print(HCLEN_ORDENADO)





            # -----------------------------------------------------------------
            # SEMANA 2, PONTO 3
            # TRANSFORMAR OS COMPRIMENTOS EM CODIGOS DE HUFFMAN DECIMAIS
            # -----------------------------------------------------------------

            CODIGOS_BINARIO = codigosHuffman(HCLEN_ORDENADO, 19)
            # print("CODIGOS HUFFMAN (BINARIO) -> ", CODIGOS_BINARIO)





            # -----------------------------------------------------------------
            # SEMANA 3, PONTO 4
            # LE E ARMAZENA NUM ARRAY OS HLIT + 257 COMPRIMENTOS DOS CÓDIGOS
            # -----------------------------------------------------------------

            # LITERAIS / COMPRIMENTOS
            CODLITCOMP = defBitABit(self, HLIT, CODIGOS_BINARIO)





            # -----------------------------------------------------------------
            # SEMANA 4, PONTO 5
            # LE E ARMAZENA NUM ARRAY OS HDIST + 1 COMPRIMENTOS DE CODIGO
            # -----------------------------------------------------------------

            # COMPRIMENTOS / DISTANCIAS
            CODDIST = defBitABit(self, HDIST, CODIGOS_BINARIO)





            # -----------------------------------------------------------------
            # SEMANA 5, PONTO 6
            # CODIGOS DE HUFFMAN DOS DOIS ALFABETOS (literais / comprimentos e distancias)
            # -----------------------------------------------------------------

            CODIGOS_CODLITCOMP_BINARIO = codigosHuffman(CODLITCOMP, HLIT)
            CODIGOS_CODDIST_BINARIO = codigosHuffman(CODDIST, HDIST)

            # print(CODIGOS_CODDIST_BINARIO)
            # print(CODIGOS_CODLITCOMP_BINARIO)





            # -----------------------------------------------------------------
            # SEMANA 5, PONTO 7
            # DESCODIFICACAO DE DADOS ATRAVES DAS DUAS ARVORES
            # -----------------------------------------------------------------

            # agoritmo para o deflate (Doc1 - 46/49)
            OUTPUT_STREAM = []

            simbolo = 0
            strLer = ""

            # enquanto o simbolo não é o fim do bloco
            while simbolo != 256:
                # lê bit a bit até encontra""r um simbolo LIT
                strLer = strLer + str(self.readBits(1))

                if strLer in CODIGOS_CODLITCOMP_BINARIO:
                    # simbolo LIT encontrado!
                    simbolo = CODIGOS_CODLITCOMP_BINARIO.index(strLer)

                    # se simbolo < 256 -> literal -> adicionar diretamenta no OutputStream
                    if simbolo < 256:
                        # print("pos encontrado, simbolo ", simbolo)
                        OUTPUT_STREAM.append(simbolo)
                        strLer = ""
                    else:
                        # se nao

                        # fim do bloco encontrado -> acabar ciclo while
                        if simbolo == 256:
                            break
                        else:
                            # simbolo >256 encontrado! -> simbolo especial -> tratamento
                            # print('pos especial, simbolo ', simbolo)

                            strLer = ""

                            # tratamento especial recorrendo ás tabelas Doc1 e Arrays auxiliares criados acima
                            if simbolo < 265:
                                comp = (simbolo - 257) + 3

                            else:
                                dif = simbolo - 265
                                bitsALer = BITSHLIT[dif]
                                lenComp = LENGTHSHLIT[dif]

                                comp = lenComp + self.readBits(bitsALer)

                            # encontrar dist recorrendo aos arrays auxiliares (BITHDIST e LENGHT)
                            strLerDist = str(self.readBits(1))
                            while strLerDist not in CODIGOS_CODDIST_BINARIO:
                                strLerDist = strLerDist + str(self.readBits(1))

                            distEspecial = CODIGOS_CODDIST_BINARIO.index(strLerDist)

                            bitsALer = BITSHDIST[distEspecial]
                            lengthPercorrer = LENGTHSHDIST[distEspecial]

                            dist = lengthPercorrer + self.readBits(bitsALer)

                            # depois de dist e length encontrados calcular o simbolo definitivo a inserir no outputStream (<length,dist>)
                            for i in range(comp):
                                OUTPUT_STREAM.append(OUTPUT_STREAM[len(OUTPUT_STREAM) - dist])

                            # print('-> length', comp)
                            # print('-> dist', distEspecial)





            # -----------------------------------------------------------------
            # SEMANA 6, PONTO 8
            # ESCREVER FICHEIRO
            # -----------------------------------------------------------------

            # abrir o ficheiro em modo de escrita binária
            f = open(self.gzh.fName, 'wb')

            # escrever os bytes do outputstream para o ficheiro
            f.write(bytes(OUTPUT_STREAM))

            # fechar o ficheiro
            f.close()

            # update number of blocks read
            numBlocks += 1

        # close file

        self.f.close()
        print("End: %d block(s) analyzed." % numBlocks)

    def getOrigFileSize(self):
        ''' reads file size of original file (before compression) - ISIZE '''

        # saves current position of file pointer
        fp = self.f.tell()

        # jumps to end-4 position
        self.f.seek(self.fileSize - 4)

        # reads the last 4 bytes (LITTLE ENDIAN)
        sz = 0
        for i in range(4):
            sz += self.f.read(1)[0] << (8 * i)

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

        mask = (2 ** n) - 1
        value = self.bits_buffer & mask

        if not keep:
            self.bits_buffer >>= n
            self.available_bits -= n

        return value


if __name__ == '__main__':

    # gets filename from command line if provided
    filename = "FAQ.txt.gz"
    if len(sys.argv) > 1:
        fileName = sys.argv[1]

    # decompress file
    gz = GZIP(filename)
    gz.decompress()
