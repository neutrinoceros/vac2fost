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
        res = func(*args, **kwargs)
        print("ok")
        return res
    return modfunc

@mydecok("function 1")
def f():
    sleep(5)
@mydecok("function 2")
def ff(a:int):
    sleep(a)
@mydecok("function 3")
def fff(a, b=6):
    sleep(5)
    return a**b
# "script"
f()
ff(2)
res = fff(5)

assert res == 15625
