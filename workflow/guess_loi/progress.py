from os.path import devnull


class Progress:
    def __init__(self, filename):
        if filename is None:
            filename = devnull
        self._filename = filename

    def __enter__(self):
        self._out = open(self._filename, "wt")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._out.close()

    def track(self, label, items, n=None):
        if n is None:
            items = list(items)
            n = len(items)

        if n == 0:
            return

        self._write('T%s\n' % label)

        for i, item in enumerate(items):
            self._out.write('P%d\n' % round(i/n * 100))
            self._out.flush()
            yield item

        self.complete()

    def _write(self, message):
        self._out.write(message)
        self._out.flush()

    def start(self, label):
        self._write('U%s\n' % label)

    def complete(self):
        self._write('P100\n')
