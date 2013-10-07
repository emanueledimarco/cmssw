
class change:
    def __init__(self,ch):
        self.baseW = ctypes.c_float(0.0565850053)
        self.b2 = ch.GetBranch("baseW")
        self.b2.SetAddress(self.baseW)

    def change(self,ch):
#         self.baseW.value = random.random()
        return True


# EOF
