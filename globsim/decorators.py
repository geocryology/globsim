def check(f):
    def wrapper(*args, **kwargs):
        if args[0].skip_checks:
            return
        else:
            return f(*args, **kwargs)
    
    return wrapper