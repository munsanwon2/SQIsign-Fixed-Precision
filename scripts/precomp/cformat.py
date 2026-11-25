#!/usr/bin/env python3
from sage.all import Integer

from parameters import ibz_limbs, p_bits, ref_radix


class Ibz:
    def __init__(self, value):
        self.value = int(value)

    def limbs(self):
        limit = 1 << (64 * ibz_limbs - 1)
        if not -limit <= self.value < limit:
            raise OverflowError(f"{self.value} does not fit a signed {64 * ibz_limbs}-bit ibz_t")
        mask = (1 << 64) - 1
        return [(self.value >> (64 * i)) & mask for i in range(ibz_limbs)]

    def _literal(self):
        literals = []
        for word in self.limbs():
            if word == 0:
                literals.append("0")
            elif word == (1 << 64) - 1:
                literals.append("UINT64_MAX")
            else:
                literals.append(f"UINT64_C({word:#x})")
        rows = [", ".join(literals[i:i + 12]) for i in range(0, len(literals), 12)]
        return "{\n    " + ",\n    ".join(rows) + "\n}"


class FpEl:
    def __init__(self, n, p, montgomery=True):
        self.n = n
        self.p = p
        self.montgomery = montgomery

    def _literal(self, word_size, broadwell=False):
        radix = word_size if broadwell else ref_radix[word_size]
        limbs = (p_bits + radix - 1) // radix
        factor = 1 << (radix * limbs) if self.montgomery else 1
        el = int((self.n * factor) % self.p)
        words = [(el >> (radix * i)) & ((1 << radix) - 1) for i in range(limbs)]
        literal = "{" + ", ".join(map(hex, words)) + "}"
        return "{" + literal + "}" if broadwell else literal


class Object:
    def __init__(self, ty, name, obj):
        if '[' in ty:
            idx = ty.index('[')
            depth = ty.count('[]')

            def shape(values, depth):
                if not depth:
                    return ()
                assert isinstance(values, (list, tuple))
                inner, = {shape(value, depth - 1) for value in values}
                return (len(values),) + inner

            dims = shape(obj, depth)
            self.ty = ty[:idx], ''.join(f'[{d}]' for d in dims)
        else:
            self.ty = ty, ''
        self.name = name
        self.obj = obj

    def _declaration(self):
        return f'extern const {self.ty[0]} {self.name}{self.ty[1]};'

    def _literal(self):
        def rec(obj):
            if isinstance(obj, (int, Integer)):
                return str(obj) if obj < 256 else hex(obj)
            if isinstance(obj, Ibz):
                return obj._literal()
            if isinstance(obj, FpEl):
                return (
                    "\n#if defined(SQISIGN_GF_IMPL_BROADWELL)\n"
                    + obj._literal(64, broadwell=True)
                    + "\n#elif RADIX == 64\n" + obj._literal(64)
                    + "\n#else\n" + obj._literal(32)
                    + "\n#endif\n"
                )
            if isinstance(obj, (list, tuple)):
                return '{' + ', '.join(map(rec, obj)) + '}'
            if isinstance(obj, str):
                return obj
            raise TypeError(f'unknown type {type(obj)} in Formatter')

        return rec(self.obj)

    def _definition(self):
        return f'const {self.ty[0]} {self.name}{self.ty[1]} = ' + self._literal() + ';'


class ObjectFormatter:
    def __init__(self, objs):
        self.objs = objs

    def header(self, file=None):
        for obj in self.objs:
            print(obj._declaration(), file=file)

    def implementation(self, file=None):
        for obj in self.objs:
            definition = obj._definition()
            print("\n".join(line.rstrip() for line in definition.splitlines()), file=file)
