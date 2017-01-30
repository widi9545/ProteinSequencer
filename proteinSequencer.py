import math, getopt, sys, os.path, time


def main(argv):

	#here is the main runner sequence 
	inputFile= ''
	inputSequence = ''
	inputKmer = ''
	
	try:
		opts, args = getopt.getopt(argv, "f:s:k:")

	except getopt.GetoptError:
		print "error"
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-f':
			if os.path.isfile(arg):
				inputFile = arg
			else:
				sys.exit(2)
		if opt == '-k':
			inputKmer = arg
		if opt == '-s':
			inputSequence = arg

	print "File:", inputFile
	idFinder(inputFile, inputSequence, inputKmer)



def gcCounter(sequenceCounter, sequences):
	#here is how we calculate the GC of a given sequence
	#It works simple - it matchs each character to see which of the four bases it is
	#we have a total basecount called baseCount and a GC counter called GC
	#once we encounter the next sequence marked with '>' we break out of the sequence and return
	GC = 0
	baseCount = 0

	# we take the sequenceCounter from what calculated in idFinder so that we can start at the appropriate sequence
	for x in range(sequenceCounter, len(sequences)):

		if sequences[x] == "G" or sequences[x] == "C":
			GC = GC + 1
			baseCount = baseCount + 1
		if sequences[x] == "A" or sequences [x] == "T":
			baseCount = baseCount + 1
		#we stop at the next sequence
		if sequences[x] == ">":
			break
	gcPercentage = (float(GC)/float(baseCount)) * 100
	return gcPercentage

def kmerCounter(sequenceCounter, sequences, kmer,sequenceID):
	#here is where things got a little hairy - I tried implementing the naive string pattern matcher from CLRS
	#but it didn't quiiiite work, and i'm a little off on my numbers. It works okay on the sample_5mer.fasta file that was given
	#but I dont know how well it works on the SCER_r64.fasta file that was given - especially considering the file is MASSIVE
	#i'm definitely going to be off on my numbers for this one, looking forward to seeing a solution!
	kmerCount = 0
	kmerID = 0
	currentBase = 0
	#where the sequence starts and ends for a specific kmer
	start = []
	end = []

	kmerList = []
	
	#again, we take the sequence counter from the idFinder so that we can start at the appropriate sequence
	for x in range(sequenceCounter, len(sequences)):
		#we need to match on the pattern of the kmer 
		if sequences[x] == kmer[kmerID]:
			#if there is a match on a character, we increase the "current base" counter

			currentBase = currentBase + 1
			if kmerID < (len(kmer)-1):
				kmerID = kmerID + 1
		#if the "CurrentBase" counter matches the length of the kmer, we increase the kmerCount and set the currentBase counter back to 0 - this is primarily where my bug comes into play, but it seems to be minimal and only in certain sequences
		if currentBase == len(kmer):
			start.append(x-len(kmer))
			end.append(x)
			currentBase = 0
			kmerCount = kmerCount + 1
		else:
			kmerID = 0
		#again we break at the next sequence
		if sequences[x] == ">":
			break
	for i in range(0, len(start)):
		print sequenceID, "Diment", "match", start[i], end[i], "100.", "-", ".", "ID =",i



def idFinder(fileToBeRead, sequenceID, kmer):
	#here is where we actually find the sequence!

	fo = open(fileToBeRead, 'r')
	sequences = []
	output = []
	sequenceCounter = 1
	idCounter = 0
	#when we wish to find a specific kmer, we store it here
	kmerInput = kmer

	fileOpen = fo.read()

	#we append the entirety of the file to a list
	for x in fileOpen:
		sequences.append(x)

	#here we look for the specific sequence ID that we want
	#we append each character of the sequence ID to a list we call output - where we output the ID
	#once the output matches the length of the sequence ID, we return
	#all the while, we increment a "Sequence Counter" to keep track of the index of the specific sequence we want
	for x in range(0, len(sequences)):
		if sequences[x] == sequenceID[idCounter]:
			if idCounter < (len(sequenceID) - 1):
				idCounter = idCounter + 1
			output.append(sequences[x])
			if len(output) == len(sequenceID):
				break
		else:
			idCounter = 0
			output = []
		
		sequenceCounter = x + 1
	#we calculate the GC here
	GC = gcCounter(sequenceCounter, sequences)
	#here is where we locate specific kmers
	print "Seq: ", sequenceID
	print "K-Mer: ", kmer
	numKmer = kmerCounter(sequenceCounter, sequences, kmerInput,sequenceID)
	

	
if __name__ == "__main__":
	main(sys.argv[1:])
	
