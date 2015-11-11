
class R_environment(object):
    def __init__(self):
        import rpy2.robjects, rpy2.robjects.numpy2ri
        rpy2.robjects.numpy2ri.activate()
        self.r = rpy2.robjects.r
        self.r("library(vcd); library(ggplot2)")
        
        self.env = self.r("new.env()")
    
    def __call__(self, text):
        return self.r.eval(self.r.parse(text=text),envir=self.env)
    
    def __setitem__(self, key, val):
        self.env[key] = val
    
    def __getitem__(self, key):
        return self.env[key]
    
    def the(self, text):
        result = self(text)
        assert len(result) == 1, "Expected 1 value but got %d in %s" % (len(result), text)
        return result[0]

    #def get_p(self, expr, term):
    #    self["x_expr"] = self(expr)
    #    self["x_term"] = term
    #    if term in self("rownames(x_expr$coefficients)"):
    #        return self("x_expr$coefficients[x_term,'Pr(>|z|)']")[0]
    #    else:
    #        return 1.0


if __name__ == "__main__":
    import numpy
    env = R_environment()
    env["x"] = numpy.array([1,2,3])
    env("print(x)")