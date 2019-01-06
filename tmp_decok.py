from time import sleep

def parameterized(dec):
    """source: https://stackoverflow.com/questions/5929107/decorators-with-parameters"""
    def layer(*args, **kwargs):
        def repl(f):
            return dec(f, *args, **kwargs)
        return repl
    return layer

@parameterized
def mydecok(func, mess):
    def modfunc(*args, **kwargs):
        print(mess.ljust(50), end="... ", flush=True)
        func(*args, **kwargs)
        print("ok")
    return modfunc

@mydecok("hello")
def fff():
    sleep(5)

fff()
