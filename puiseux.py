r"""
Puiseux series expansions by dynamic evaluation (experimental)


Authors: Ariane Carrance, Marc Mezzarobba, 2023--2024

Inspired in part by gfun:-algeqtoseries by Bruno Salvy
(https://perso.ens-lyon.fr/bruno.salvy/software/the-gfun-package/).


EXAMPLES::

    sage: from puiseux import puiseux
    sage: P.<x> = QQ[]; Q.<y> = P[]

Duval 1989::

    sage: F = (x^2+y^2)^3 - 4*x^2*y^2
    sage: f = puiseux(F, 4)
    sage: f[0]
    -1/2*x^2 - 3/16*x^4 - 39/256*x^6 - 323/2048*x^8 + O(x^10)
    sage: f[1]
    1/2*x^2 + 3/16*x^4 + 39/256*x^6 + 323/2048*x^8 + O(x^10)

    sage: F = y^16-4*y^12*x^6-4*y^11*x^8+y^10*x^10+6*y^8*x^12+8*y^7*x^14+14*y^6*x^16+4*y^5*x^18+y^4*(x^20-4*x^18)-4*y^3*x^20+y^2*x^22+x^24  # p. 140
    sage: f = puiseux(F, 2); f
    [-x^(3/2) + alg1*x^(7/4) + O(x^(9/4)),
    alg01*x^(3/2) + alg2*x^(7/4) + O(x^(9/4))]
    sage: f[0].base_ring()
    Univariate Quotient Polynomial Ring in alg1 over Rational Field with modulus 16*alg1^4 + 4*alg1^2 + 1
    sage: f[1].base_ring()
    Univariate Quotient Polynomial Ring in alg2 over Univariate Quotient Polynomial Ring in alg01 over Rational Field with modulus alg01^3 - alg01^2 + alg01 - 1 with modulus 256*alg2^4 - 64*alg01*alg2^2 + 16*alg01^2

    sage: G = y^16-2*y^13*x^5-4*y^12*x^6+y^10*x^10+2*y^9*x^11+y^8*(8*x^13+6*x^12)-2*y^6*x^16+y^5*(-4*x^18+2*x^17)+y^4*(-x^20+8*x^19-4*x^18)+y^2*x^22-2*y*x^23+x^2  # p. 144 (not checked)
    sage: f = puiseux(G, 2)
    sage: f[0]
    alg0*x^(1/8) - 1/8*alg0^14*x^(19/4) + O(x^(45/8))
    sage: f[0].base_ring()
    Univariate Quotient Polynomial Ring in alg0 over Rational Field with modulus alg0^16 + 1

Here one of the coefficients of the expansion vanishes for some of the possible
values of ``alg0``::

    sage: from puiseux import puiseux
    sage: P.<x> = QQ[]; Q.<y> = P[]
    sage: p = ((y^2 - 2*x + 3*x^5)*(y^2 + 2*x + 3*x^3))
    sage: f = puiseux(p, 0); f
    [O(x^(1/2)), O(x^(1/2)), O(x^(1/2)), O(x^(1/2))]
    sage: f[0].base_ring()
    Rational Field
    sage: f = puiseux(p, 1); f
    [alg00*x^(1/2) + O(x^(9/2)), alg01*x^(1/2) + O(x^(5/2))]
    sage: f[0].base_ring()
    Univariate Quotient Polynomial Ring in alg00 over Rational Field with modulus alg00^2 - 2
    sage: f[1].base_ring()
    Univariate Quotient Polynomial Ring in alg01 over Rational Field with modulus alg01^2 + 2
    sage: puiseux(p, 4)
    [alg00*x^(1/2) - 3/4*alg00*x^(9/2) - 9/32*alg00*x^(17/2) - 27/128*alg00*x^(25/2) + O(x^(33/2)),
    alg01*x^(1/2) + 3/4*alg01*x^(5/2) - 9/32*alg01*x^(9/2) + 27/128*alg01*x^(13/2) + O(x^(17/2))]

An example from Bouttier & Carrance (EJC, 2021,
<https://www.combinatorics.org/ojs/index.php/eljc/article/view/v28i3p21>)
where the base ring contains a parameter::

    sage: Pol.<z, t, B, V, W> = QQ[]
    sage: sys = [V*(2*B*V^2*W - V^2*W^2 - B*V*W - V^2*W + 2*B*V - 2*V*W - B - 2*V + W + 1),
    ....:        2*V^2 - V + t,
    ....:        V^2*W^2 + 2*V*W + W*z - W + z]
    sage: Id = ideal(sys)
    sage: El = Id.elimination_ideal([V,W])
    sage: pol = El.gen(0)//(t*(1+t))
    sage: trans = pol(t=(1-t)/8).numerator()
    sage: trans = Frac(QQ['z'])['t']['B'](trans)
    sage: trans
    ((-16*z^2 + 32*z - 16)*t^3 + (48*z^2 - 96*z + 48)*t^2 + (-48*z^2 +
    96*z - 48)*t + 16*z^2 - 32*z + 16)*B^4 + ((32*z^2 - 64*z + 32)*t^3 +
    (-256*z^3 + 416*z^2 - 64*z - 96)*t^2 + (512*z^3 - 928*z^2 + 320*z
    + 96)*t - 256*z^3 + 480*z^2 - 192*z - 32)*B^3 + ((-24*z^2 + 40*z -
    16)*t^3 + (256*z^3 - 600*z^2 + 296*z + 48)*t^2 + (-1024*z^4 + 1536*z^3
    + 312*z^2 - 776*z - 48)*t + 1024*z^4 - 1792*z^3 + 312*z^2 + 440*z +
    16)*B^2 + ((8*z^2 - 8*z)*t^3 + (-64*z^3 + 200*z^2 - 136*z)*t^2 +
    (-1152*z^3 + 792*z^2 + 360*z)*t + 1728*z^3 - 1512*z^2 - 216*z)*B -
    z^2*t^3 + 27*z^2*t^2 - 243*z^2*t + 729*z^2
    sage: [f] = puiseux(trans, 4); f
    alg0 + (((z^2 - 5/8*z)/(z^2 - 1/2*z + 1/16))*alg0 + 9/16*z/(z^2 - 1/2*z +
    1/16))*t + alg1*t^(3/2) + (((z^4 - 9/8*z^3 + 47/128*z^2 - 7/128*z)/(z^4 -
    z^3 + 3/8*z^2 - 1/16*z + 1/256))*alg0 + (19/32*z^3 - 55/256*z^2 +
    11/256*z)/(z^4 - z^3 + 3/8*z^2 - 1/16*z + 1/256))*t^2 + O(t^(5/2))
    sage: f.base_ring()
    Univariate Quotient Polynomial Ring in alg1 over Univariate Quotient
    Polynomial Ring in alg0 over Fraction Field of Univariate Polynomial
    Ring in z over Rational Field with modulus (z - 1)*alg0^2 + (-8*z^2 +
    7*z + 1)*alg0 - 27/4*z with modulus (1024*z^4 - 1792*z^3 + 960*z^2 -
    208*z + 16)*alg1^2 - 64*z^2
"""


import itertools

from sage.arith.misc import gcd
from sage.categories.fields import Fields
from sage.geometry.newton_polygon import NewtonPolygon
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.rings.polynomial.polynomial_quotient_ring import PolynomialQuotientRing_generic
from sage.rings.polynomial.polynomial_quotient_ring_element import PolynomialQuotientRingElement
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.puiseux_series_ring import PuiseuxSeriesRing
from sage.rings.quotient_ring import QuotientRing_generic
from sage.structure.richcmp import op_EQ, op_NE
from sage.structure.unique_representation import UniqueRepresentation


class DynamicExtensionSplit(ZeroDivisionError):

    def __init__(self, proj):
        assert len(proj) >= 2
        self.domain = proj[0].domain()
        self.proj = proj
        super().__init__(f"Reducible extension: {self.domain}")
        assert all(mor.domain() is self.domain for mor in proj)


class DynamicExtensionElement(PolynomialQuotientRingElement):

    # TODO: conversion to plain PolynomialQuotientRingElement's

    def __invert__(self):
        if self._polynomial.is_zero():
            raise ZeroDivisionError()
        if self._polynomial.is_one():
            return self
        parent = self.parent()
        if self._polynomial.is_unit():
            inv_pol = self._polynomial.inverse_of_unit()
            return self.__class__(parent, inv_pol)
        modulus = parent.modulus()
        g, _, a = modulus.xgcd(self._polynomial)
        if g.degree() == 0:
            return self.__class__(self.parent(), (~g[0])*a, check=False)
        else:
            # We discovered a proper factor g of modulus
            raise DynamicExtensionSplit(parent.split(g))

    def is_zero(self):
        if self._polynomial.is_zero():
            return True
        elif self._polynomial.degree() == 0:
            return self._polynomial[0].is_zero()
        g = self.parent().modulus().gcd(self._polynomial)
        if g.degree() == 0:
            # if we are working over another dynamic extension, just computing
            # the degree should have been enough to split it if needed
            assert g[0].is_unit()
            return False
        else:
            raise DynamicExtensionSplit(self.parent().split(g))

    def _richcmp_(self, other, op):
        eq = (self - other).is_zero()
        if op == op_EQ:
            return eq
        elif op == op_NE:
            return not eq
        else:
            raise ValueError("no order is defined on {self.parent()}")

    def is_unit(self):
        return self != 0

    def __int__(self):
        raise NotImplementedError


class DynamicExtension(UniqueRepresentation, PolynomialQuotientRing_generic):
    r"""
    TESTS::

        sage: from puiseux import DynamicExtension, DynamicExtensionSplit
        sage: Pol.<x> = QQ[]
        sage: E = DynamicExtension(Pol, (x^2-2)*(x^2-3), 'a')
        sage: a = E.gen()
        sage: a^2-2
        a^2 - 2
        sage: a^2 == 0
        False
        sage: a^2 != 0
        True
        sage: ~a*a
        1
        sage: a^2 == 2
        Traceback (most recent call last):
        ...
        DynamicExtensionSplit: Reducible extension: Univariate Quotient
        Polynomial Ring in a over Rational Field with modulus x^4 - 5*x^2 + 6
        sage: try:
        ....:     (a^3 - 2*a).is_zero()
        ....: except DynamicExtensionSplit as exn:
        ....:     tuple(mor(a^3 - 2*a) for mor in exn.proj)
        (0, a1)
    """

    Element = DynamicExtensionElement

    def __init__(self, ring, polynomial, name=None, category=None):
        self._PolynomialQuotientRing_generic__ring = ring
        self._PolynomialQuotientRing_generic__polynomial = polynomial

        if not polynomial.leading_coefficient().is_unit():
            raise ValueError("polynomial must have unit leading coefficient")

        QuotientRing_generic.__init__(self, ring, ring.ideal(polynomial), names=name, category=Fields())
        self._base = ring  # backwards compatibility -- different from QuotientRing_generic

    def is_field(self, proof=True):
        return True

    def split(self, fac0):
        fac1, rem = self.modulus().quo_rem(fac0)
        if not rem.is_zero():
            raise ValueError(f"{fac0} is not a factor of {self.modulus()}")
        def mor(i, fac):
            if fac.degree() == 1:
                rt = -fac[0]/fac[1]
            else:
                name = self.variable_name() + str(i)
                fac = fac.change_variable_name(name)
                rt = self.__class__(fac.parent(), fac, name).gen()
            return self.hom([rt], codomain=rt.parent())
        return (mor(0, fac0), mor(1, fac1))


# would it be useful to accept a minimum valuation instead of the boolean flag
# "positive_valuation"?
def puiseux(pol, order, positive_valuation=False, *, used_names=None, depth=0,
            verbose=False):
    r"""
    Compute a complete set of Puiseux series expansions of solutions of pol at
    the origin.

    INPUT:

    - ``pol`` - TODO
    - ``order`` - number of terms

    OUTPUT:

    A list of Puiseux series, each with coefficients in a (potentially
    different) extension of the base ring of ``pol``, with the property that
    considering all embeddings of the coefficient rings in an algebraic closure
    of the base ring yields a full set of solutions over that algebraic
    closure.
    """

    def dbg(msg):
        if verbose:
            print("  "*depth + msg)

    Pol_xy = pol.parent()
    y = Pol_xy.gen()
    Pol_x = Pol_xy.base_ring()
    Csts = Pol_x.base_ring()
    Puiseux = PuiseuxSeriesRing(Csts, Pol_x.variable_name())

    Pol_y = PolynomialRing(Csts, y)

    if used_names is None:
        used_names = set(Pol_xy.variable_names_recursive())

    # Exact solutions
    valuation = pol.valuation()  # in y
    sol = [Puiseux.zero()]*valuation
    pol >>= valuation

    polygon = NewtonPolygon([(i, c.valuation()) for i, c in enumerate(pol)])
    slopes = polygon.slopes(repetition=False)
    # in recursive calls, nonpositive valuations correspond to known terms
    if positive_valuation:
        slopes = [slope for slope in slopes if slope < 0]
    for slope in slopes:
        pol1 = pol
        p, q = -slope.numerator(), slope.denominator()
        # Reduce to a horizontal edge:
        # x[new] = x[old]^(1/q)
        #   (--> y[old] = cst·x[old]^(p/q) + ··· = x[new]^p + ···)
        # y[new] = x[old]^(-p/q)·y[old] = x[new]^(-p)·y[old]
        x = Pol_x.gen()
        xx = LaurentPolynomialRing(Csts, Pol_x.variable_name()).gen()
        pol1 = Pol_xy([c(x**q) for c in pol1])(xx**p*y)
        minval = min(c.valuation() for c in pol1)
        pol1 = Pol_xy(xx**(-minval)*pol1)
        eq_cst_coeff = Pol_y([c[0] for c in pol1])
        dbg(f"{slope=}")# {pol1=} {eq_cst_coeff=}")

        # Ignore zero roots since we are looking for a solution of valuation
        # exactly p/q.  This yields a polynomial of degree equal to the length
        # of the current edge of the Newton polygon.
        eq_cst_coeff >>= eq_cst_coeff.valuation()

        if order == 0:
            # Return len(slope) O(·) terms defined over Csts. We could also
            # continue a bit further and return O(·) terms defined over each
            # Ext that stand for deg(Ext) different tails starting with a
            # constant belonging to Ext, but this seems to make the code more
            # complicated for little benefit.
            # XXX Should come before the computation of eq_cst_coeff and use
            # the length of the slope instead.
            big_oh_term = big_oh(Puiseux, p, q)
            sol.extend([big_oh_term]*eq_cst_coeff.degree())
            continue

        decomp = eq_cst_coeff.squarefree_decomposition()
        for sqf, mult in decomp:
            sqf = _remove_content(sqf)
            if sqf.degree() == 1:
                cst_terms = [-sqf[0]/sqf[1]]
                # Force a split if possible
                # cst_terms[0].is_zero()
            else:
                alg = _choose_name(used_names, 'alg')
                Pol_alg = Pol_y.change_var(alg)
                cst_terms = [DynamicExtension(Pol_alg, Pol_alg(sqf), alg).gen()]
            while cst_terms:
                y0 = cst_terms.pop()
                Ext = y0.parent()
                dbg(f"{Ext=}")
                pol1_Ext = change_constants(pol1, Ext)
                y_Ext = change_constants(y, Ext)
                try:
                    tails = puiseux(pol1_Ext(y0 + y_Ext), order - 1,
                                    positive_valuation=True,
                                    used_names=used_names,
                                    depth=depth + 1,
                                    verbose=verbose)
                except DynamicExtensionSplit as split:
                    if split.domain is Ext:
                        # XXX Instead of redoing the whole computation, we
                        # could maybe apply the morphisms to the terms computed
                        # up to this point... (But how exactly do we handle
                        # recursive splits and the like then?)
                        dbg("splitting Ext")
                        cst_terms.extend(proj(y0) for proj in reversed(split.proj))
                    else:
                        raise
                else:
                    # Note that this coerces y0 into PuiseuxSeries(Ext).
                    sol.extend(ramify((y0 + y1).shift(p), q) for y1 in tails)
    return sol

def change_constants(pol, Csts):
    Pol_xy = pol.parent()
    Pol_x = pol.base_ring()
    if Pol_x.base_ring() is Csts:
        return pol
    NewPol_x = Pol_x.change_ring(Csts)
    NewPol_xy = Pol_xy.change_ring(NewPol_x)
    new_pol = NewPol_xy([NewPol_x(c) for c in pol])
    return new_pol

def _choose_name(names, stem):
    for i in itertools.count():
        candidate = stem + str(i)
        if candidate in names:
            continue
        names.add(candidate)
        return candidate

def ramify(f, index):
    Puiseux = f.parent()
    return Puiseux(f.laurent_part(), e=index*f.ramification_index())

def big_oh(Puiseux, p, q):
    return Puiseux(Puiseux.laurent_series_ring().zero().add_bigoh(p), e=q)

def _remove_content(pol):
    Base = pol.base_ring()
    pol = pol.numerator()
    # TODO: also remove integer content in 2 variables, etc.
    if isinstance(Base, DynamicExtension):
        g = gcd(c._polynomial for c in pol)
        cst = Base.element_class(Base, g)
    else:
        # over QQ, this works, while content() returns 1
        cst = gcd(c for c in pol)
    return ~cst*pol
