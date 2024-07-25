class LinkedLinePDE:
    def __init__(self, end_d, a):
        self.end_d = end_d
        self.a = a
        self.next = None
        self.pre = None

    # 将一个方法转换为只读属性
    @property
    def start_k(self):
        return self.pre.end_d if self.pre is not None else None

    def link_forward(self, line):
        self.next = line
        line.pre = self
