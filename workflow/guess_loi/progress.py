class Progress:
    def __init__(self, filename):
        self._filename = filename

    def __enter__(self):
        self._out = open(self._filename, "wt")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._out.close()

    def step(self, label):
        self._out.write("S %s\n" % label)

    def track(self, items):
        items = list(items)
        n = len(items)
        for i, item in enumerate(items):
            self._out.write('P %d\n' % round(i/n * 100))
            self._out.flush()
            yield item

        self.complete()

    def complete(self):
        self._out.write('P 100\n')
        self._out.flush()
