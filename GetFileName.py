#class to get file name

class FileName:
	
	def __init__(self):
		self
		
	def GetFileName(self):
		
		import os, sys
		filenames = []
		for i in os.listdir(os.getcwd()):
			if i.endswith(".inp"):
				filenames.append(i)
			
			else:
				continue
				
		filenames.sort()	
		self.listfn=filenames