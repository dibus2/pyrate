try :
	import logging
except ImportError :
	exit("error while loading logging module")
#################
## Logging system
#################

def logger(func):
	def inner(mess,**keyargs):
		for key,val in keyargs.items():
			if key == 'verbose' and val:
				print(mess)
		return func(mess)
	return inner

@logger
def loggingInfo(mess):
	logging.info(mess)
@logger
def loggingDebug(mess):
	logging.debug(mess)
@logger
def loggingCritical(mess):
	logging.critical(mess)

