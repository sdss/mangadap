#Expose rundap module
try:
    from mangadap.survey.rundap import rundap
except ImportError as e:
    print(e)
    print('WARNING: rundap will be unavailable!')
