[loggers]
keys=root,sampleLogger

[handlers]
keys=consoleHandler

[formatters]
keys=sampleFormatter

[logger_root]
level=INFO
handlers=consoleHandler

[logger_sampleLogger]
level=INFO
handlers=consoleHandler
qualname=sampleLogger
propagate=0

[handler_consoleHandler]
class=StreamHandler
level=INFO
formatter=sampleFormatter
args=(sys.stdout,)

[formatter_sampleFormatter]
format=%(asctime)s.%(msecs)03d|%(levelname)7s|%(name)s|%(message)s
datefmt=%Y-%m-%dT%H:%M:%S 
