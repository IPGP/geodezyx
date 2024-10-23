# THE loggers KEY MUST DESCRIBES THE CURRENT MODULE/PROJECT!!!!

###### USEFUL TUTOS
## https://coderzcolumn.com/tutorials/python/logging-simple-guide-to-log-events-in-python
## https://coderzcolumn.com/tutorials/python/logging-config-simple-guide-to-configure-loggers-from-dictionary-and-config-files-in-python#2

log_config_dict = {
    "version":1,
    'disable_existing_loggers': False,

    
    "root":{
        "handlers" : ["console_root"],
        "level": "WARN",
        "propagate": False
    },
    
    
    'loggers': {
        "geodezyx" : {  # THIS KEY MUST DESCRIBES THE CURRENT MODULE!!!!
            "handlers" : ["console_gyxz"],
            "level": "DEBUG",
            "propagate": False

        }
    },
    
    "handlers":{        
        "console_root":{
#            "formatter":"fmtgzyx1",
            "class":"logging.StreamHandler",
            "level":"DEBUG",
        },
        "console_gyxz":{
            "formatter":"fmtgzyx1",
            "class":"logging.StreamHandler",
            "level":"DEBUG",
        }
    },
    "formatters":{
        "fmtgzyx1": {
             ##### n-sized fct name space, truncated after n (here n = 15)
            "fmt": "%(asctime)s.%(msecs)03d|%(log_color)s%(levelname).1s%(reset)s|%(log_color)s%(funcName)-15s%(reset)s|%(message)s",
            ##### exact fct name space (<= n), truncated after n (here n = 15)
            #"format": "%(asctime)s.%(msecs)03d|%(log_color)s%(levelname).1s%(reset)s|%(log_color)s%(funcName).15s%(reset)s|%(message)s",
            '()': 'colorlog.ColoredFormatter',
            "datefmt":"%y%m%dT%H:%M:%S",
            "log_colors":{
                'DEBUG':    'cyan',
                'INFO':     'green',
                'WARNING':  'yellow',
                'ERROR':    'red',
                'CRITICAL': 'red,bg_white',
                },
        }
    },
}



#### LEGACY CONFIG FILE

# [loggers]
# keys=root,sampleLogger

# [handlers]
# keys=consoleHandler

# [formatters]
# keys=sampleFormatter

# [logger_root]
# level=DEBUG
# handlers=consoleHandler

# [logger_sampleLogger]
# level=DEBUG
# handlers=consoleHandler
# qualname=sampleLogger
# propagate=0

# [handler_consoleHandler]
# class=StreamHandler
# level=DEBUG
# formatter=sampleFormatter
# args=(sys.stdout,)

# [formatter_sampleFormatter]
# format=%(asctime)s.%(msecs)03d|%(levelname).1s|%(name)s|%(message)s
# datefmt=%y%m%dT%H:%M:%S 



