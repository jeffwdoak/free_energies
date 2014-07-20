#!/usr/bin/python

class testes(dict):
    def __init__(self):
        #self.update(a=1)
        #self.update(b=2)
        self.a=1
        self.b=2

    def __getattr__(self,attr):
        return self[attr]

    def __setattr__(self,attr,val):
        self[attr] = val

    def function(self):
        pass

a = testes()
print a.__dict__
print a['a']
print a['b']
