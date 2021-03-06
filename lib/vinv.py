from numpy import *


def JacInv(th0, th1, th2, d1, a1, a2):
    r1 = array([(conjugate(cos(th0)) * cos(th0) * (
        conjugate(a2 * cos(th1 + th2)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(
            th1) - conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(
            th1) + (conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
            conjugate(a2 * cos(th0) * sin(th1 + th2)) - conjugate(
                cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2)))) + conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (conjugate(a2 * sin(th0) * sin(th1 + th2)) - conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))))) * sin(th0) * sin(th1)) + sin(th0) * (
                 conjugate(sin(th0)) * (
                     conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                         conjugate(a2 * cos(th0) * sin(th1 + th2)) - conjugate(
                             cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2)))) + conjugate(
                         (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                         conjugate(a2 * sin(th0) * sin(th1 + th2)) - conjugate(
                             sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))))) * sin(
                     th0) * sin(th1) - conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * (
                 conjugate(sin(th0)) * conjugate(
                     (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) + a2 * (conjugate(
                     cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                     a2 * cos(th0) * sin(th1 + th2)) + conjugate(
                     (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                     a2 * sin(th0) * sin(th1 + th2))) * sin(
                     th2)) + conjugate(a2 * cos(th1 + th2)) * (conjugate(sin(th0)) * conjugate(
                     (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) + a2 * (conjugate(
                     cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                     cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) + conjugate(
                     (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                     sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2)))) * sin(th2)))) / (
                    conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * cos(th1) * sin(
                        th0) + a1 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * cos(th0) * cos(
                        th1) ** 2 * sin(th0) + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * cos(th0) * cos(
                        th1) * cos(th1 + th2) * sin(th0) + a1 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) ** 2 * sin(
                        th0) ** 2 + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) * cos(th1 + th2) * sin(
                        th0) ** 2 - a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * cos(th1 + th2) * sin(th1) - a1 * a2 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
                        th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
                        th1) - a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(
                        th1) - a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(th1) - conjugate(
                        sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * sin(th0) * sin(
                        th1) + conjugate(sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(
                        th0) * sin(th0) * sin(th1) - a1 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
                        th1) - a1 * conjugate(
                        sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
                        th1) + a1 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
                        th1) + a1 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
                        th1) - a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) * sin(th1) - a2 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
                        th1) - a2 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
                        th1) + a2 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(
                        th1 + th2) * sin(th0) * sin(th1) + a2 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1 + th2) * sin(
                        th0) * sin(th1) - conjugate(sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * sin(
                        th0) ** 2 * sin(th1) + conjugate(sin(th0)) * conjugate(
                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) ** 2 * sin(
                        th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(
                        th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(
                        th1) - a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th1 + th2) ** 2 * sin(th0) ** 2 * sin(
                        th1) - a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) ** 2 * sin(th0) ** 2 * sin(
                        th1) - a1 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * sin(th0) ** 3 * sin(th1) - a1 * conjugate(
                        sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * sin(th0) ** 3 * sin(th1) + a1 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(
                        th1) + a1 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(
                        th1) - a2 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(
                        th1) - a2 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) + a2 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(
                        th1) + a2 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(
                        th1) + conjugate(cos(th0)) * cos(th0) * (-(conjugate(a2 * cos(th1 + th2)) * cos(th1) * (
                        1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                            a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate(
                            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                            a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0))) + conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * cos(th1) * (1 + conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (a1 * cos(th1) + a2 * cos(
                        th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                                                               a1 * cos(th1) + a2 * cos(
                                                                                   th1 + th2)) * sin(
                        th0)) + (-(conjugate(a2 * cos(th0) * sin(th1 + th2)) * (
                        cos(th0) + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                            a1 * cos(th1) + a2 * cos(th1 + th2)))) + conjugate(
                        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                     cos(th0) + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                                         a1 * cos(th1) + a2 * cos(th1 + th2))) - (
                                     conjugate(a2 * sin(th0) * sin(th1 + th2)) - conjugate(
                                         sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2)))) * (
                                     conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                         a1 * cos(th1) + a2 * cos(th1 + th2)) + sin(th0))) * sin(th1)) - conjugate(
                        a2 * cos(th1 + th2)) * (conjugate(sin(th0)) * cos(th1) * sin(th0) * (
                        1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                            a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate(
                            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                            a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + a2 * (
                                                    conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                        cos(th0) + conjugate(
                                                            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                                                            a1 * cos(th1) + a2 * cos(th1 + th2))) + conjugate(
                                                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                        conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                                            a1 * cos(th1) + a2 * cos(th1 + th2)) + sin(th0))) * sin(
                        th2)) + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * cos(th1) * sin(th1 + th2) + a1 * a2 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) ** 2 * sin(
                        th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) ** 2 * sin(
                        th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
                        th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
                        th1 + th2) + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * sin(th0) * sin(th1 + th2) + a1 * a2 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th1) ** 2 * sin(th0) ** 2 * sin(
                        th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1) ** 2 * sin(th0) ** 2 * sin(
                        th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(
                        th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(th1 + th2)), (
                    conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * cos(th1) * sin(
                        th0) - conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * sin(th0) * sin(
                        th1) - conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) * sin(th0) * sin(th1) + conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * sin(th0) * sin(th1) + conjugate(
                        sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * sin(th0) * sin(th1) + conjugate(
                        cos(th0)) * cos(th0) * (-(
                        conjugate(a2 * cos(th1 + th2)) * conjugate(
                            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(
                            th1)) + conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th1) + (
                                                    conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                                                        -conjugate(a2 * cos(th0) * sin(th1 + th2)) + conjugate(
                                                            cos(th0) * (
                                                            a1 * sin(th1) + a2 * sin(th1 + th2)))) + conjugate(
                                                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                                        -conjugate(a2 * sin(th0) * sin(th1 + th2)) + conjugate(
                                                            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))))) * cos(
                        th0) * sin(
                        th1)) + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * sin(th2) + a2 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) * sin(th2) - conjugate(a2 * cos(th1 + th2)) * (
                        a2 * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * sin(th2) + conjugate(
                            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                            conjugate(sin(th0)) * cos(th1) * sin(th0) + a2 * conjugate(
                                cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * sin(th2)))) / (
                    conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * cos(th1) * sin(
                        th0) + a1 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * cos(th0) * cos(
                        th1) ** 2 * sin(th0) + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * cos(th0) * cos(
                        th1) * cos(th1 + th2) * sin(th0) + a1 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) ** 2 * sin(
                        th0) ** 2 + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) * cos(th1 + th2) * sin(
                        th0) ** 2 - a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * cos(th1 + th2) * sin(th1) - a1 * a2 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
                        th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
                        th1) - a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(
                        th1) - a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(th1) - conjugate(
                        sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * sin(th0) * sin(
                        th1) + conjugate(sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(
                        th0) * sin(th0) * sin(th1) - a1 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
                        th1) - a1 * conjugate(
                        sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
                        th1) + a1 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
                        th1) + a1 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
                        th1) - a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) * sin(th1) - a2 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
                        th1) - a2 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
                        th1) + a2 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(
                        th1 + th2) * sin(th0) * sin(th1) + a2 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1 + th2) * sin(
                        th0) * sin(th1) - conjugate(sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * sin(
                        th0) ** 2 * sin(th1) + conjugate(sin(th0)) * conjugate(
                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) ** 2 * sin(
                        th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(
                        th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(
                        th1) - a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th1 + th2) ** 2 * sin(th0) ** 2 * sin(
                        th1) - a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) ** 2 * sin(th0) ** 2 * sin(
                        th1) - a1 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * sin(th0) ** 3 * sin(th1) - a1 * conjugate(
                        sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * sin(th0) ** 3 * sin(th1) + a1 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(
                        th1) + a1 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(
                        th1) - a2 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(
                        th1) - a2 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) + a2 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(
                        th1) + a2 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(
                        th1) + conjugate(cos(th0)) * cos(th0) * (-(conjugate(a2 * cos(th1 + th2)) * cos(th1) * (
                        1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                            a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate(
                            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                            a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0))) + conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * cos(th1) * (1 + conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (a1 * cos(th1) + a2 * cos(
                        th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                                                               a1 * cos(th1) + a2 * cos(
                                                                                   th1 + th2)) * sin(
                        th0)) + (-(conjugate(a2 * cos(th0) * sin(th1 + th2)) * (
                        cos(th0) + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                            a1 * cos(th1) + a2 * cos(th1 + th2)))) + conjugate(
                        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                     cos(th0) + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                                         a1 * cos(th1) + a2 * cos(th1 + th2))) - (
                                     conjugate(a2 * sin(th0) * sin(th1 + th2)) - conjugate(
                                         sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2)))) * (
                                     conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                         a1 * cos(th1) + a2 * cos(th1 + th2)) + sin(th0))) * sin(th1)) - conjugate(
                        a2 * cos(th1 + th2)) * (conjugate(sin(th0)) * cos(th1) * sin(th0) * (
                        1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                            a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate(
                            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                            a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + a2 * (
                                                    conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                        cos(th0) + conjugate(
                                                            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                                                            a1 * cos(th1) + a2 * cos(th1 + th2))) + conjugate(
                                                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                        conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                                            a1 * cos(th1) + a2 * cos(th1 + th2)) + sin(th0))) * sin(
                        th2)) + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * cos(th1) * sin(th1 + th2) + a1 * a2 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) ** 2 * sin(
                        th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) ** 2 * sin(
                        th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
                        th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
                        th1 + th2) + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * sin(th0) * sin(th1 + th2) + a1 * a2 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th1) ** 2 * sin(th0) ** 2 * sin(
                        th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1) ** 2 * sin(th0) ** 2 * sin(
                        th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(
                        th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(th1 + th2)), (
                    (conjugate(a2 * cos(th1 + th2)) - conjugate(a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                        -(conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th0)) + conjugate(
                            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * sin(th0)) * (
                        conjugate(cos(th0)) * cos(th0) + conjugate(sin(th0)) * sin(th0)) * sin(th1)) / (
                    conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * cos(th1) * sin(
                        th0) + a1 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * cos(th0) * cos(
                        th1) ** 2 * sin(th0) + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * cos(th0) * cos(
                        th1) * cos(th1 + th2) * sin(th0) + a1 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) ** 2 * sin(
                        th0) ** 2 + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) * cos(th1 + th2) * sin(
                        th0) ** 2 - a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * cos(th1 + th2) * sin(th1) - a1 * a2 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
                        th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
                        th1) - a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(
                        th1) - a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(th1) - conjugate(
                        sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * sin(th0) * sin(
                        th1) + conjugate(sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(
                        th0) * sin(th0) * sin(th1) - a1 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
                        th1) - a1 * conjugate(
                        sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
                        th1) + a1 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
                        th1) + a1 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
                        th1) - a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) * sin(th1) - a2 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
                        th1) - a2 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
                        th1) + a2 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(
                        th1 + th2) * sin(th0) * sin(th1) + a2 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1 + th2) * sin(
                        th0) * sin(th1) - conjugate(sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * sin(
                        th0) ** 2 * sin(th1) + conjugate(sin(th0)) * conjugate(
                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) ** 2 * sin(
                        th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(
                        th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(
                        th1) - a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th1 + th2) ** 2 * sin(th0) ** 2 * sin(
                        th1) - a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) ** 2 * sin(th0) ** 2 * sin(
                        th1) - a1 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * sin(th0) ** 3 * sin(th1) - a1 * conjugate(
                        sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * sin(th0) ** 3 * sin(th1) + a1 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(
                        th1) + a1 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(
                        th1) - a2 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(
                        th1) - a2 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) + a2 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(
                        th1) + a2 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(
                        th1) + conjugate(cos(th0)) * cos(th0) * (-(conjugate(a2 * cos(th1 + th2)) * cos(th1) * (
                        1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                            a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate(
                            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                            a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0))) + conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * cos(th1) * (1 + conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (a1 * cos(th1) + a2 * cos(
                        th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                                                               a1 * cos(th1) + a2 * cos(
                                                                                   th1 + th2)) * sin(
                        th0)) + (-(conjugate(a2 * cos(th0) * sin(th1 + th2)) * (
                        cos(th0) + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                            a1 * cos(th1) + a2 * cos(th1 + th2)))) + conjugate(
                        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                     cos(th0) + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                                         a1 * cos(th1) + a2 * cos(th1 + th2))) - (
                                     conjugate(a2 * sin(th0) * sin(th1 + th2)) - conjugate(
                                         sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2)))) * (
                                     conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                         a1 * cos(th1) + a2 * cos(th1 + th2)) + sin(th0))) * sin(th1)) - conjugate(
                        a2 * cos(th1 + th2)) * (conjugate(sin(th0)) * cos(th1) * sin(th0) * (
                        1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                            a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate(
                            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                            a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + a2 * (
                                                    conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                        cos(th0) + conjugate(
                                                            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                                                            a1 * cos(th1) + a2 * cos(th1 + th2))) + conjugate(
                                                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                        conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                                            a1 * cos(th1) + a2 * cos(th1 + th2)) + sin(th0))) * sin(
                        th2)) + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * cos(th1) * sin(th1 + th2) + a1 * a2 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) ** 2 * sin(
                        th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) ** 2 * sin(
                        th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
                        th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
                        th1 + th2) + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * sin(th0) * sin(th1 + th2) + a1 * a2 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th1) ** 2 * sin(th0) ** 2 * sin(
                        th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1) ** 2 * sin(th0) ** 2 * sin(
                        th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(
                        th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(th1 + th2)), (
                    a2 * (conjugate(a2 * cos(th1 + th2)) - conjugate(a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        sin(th0)) * (
                        -(conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th0)) + conjugate(
                            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * sin(th0)) * sin(th2)) / (
                    conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * cos(th1) * sin(
                        th0) + a1 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * cos(th0) * cos(
                        th1) ** 2 * sin(th0) + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * cos(th0) * cos(
                        th1) * cos(th1 + th2) * sin(th0) + a1 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) ** 2 * sin(
                        th0) ** 2 + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) * cos(th1 + th2) * sin(
                        th0) ** 2 - a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * cos(th1 + th2) * sin(th1) - a1 * a2 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
                        th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
                        th1) - a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(
                        th1) - a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(th1) - conjugate(
                        sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * sin(th0) * sin(
                        th1) + conjugate(sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(
                        th0) * sin(th0) * sin(th1) - a1 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
                        th1) - a1 * conjugate(
                        sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
                        th1) + a1 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
                        th1) + a1 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
                        th1) - a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) * sin(th1) - a2 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
                        th1) - a2 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
                        th1) + a2 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(
                        th1 + th2) * sin(th0) * sin(th1) + a2 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1 + th2) * sin(
                        th0) * sin(th1) - conjugate(sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * sin(
                        th0) ** 2 * sin(th1) + conjugate(sin(th0)) * conjugate(
                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) ** 2 * sin(
                        th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(
                        th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(
                        th1) - a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th1 + th2) ** 2 * sin(th0) ** 2 * sin(
                        th1) - a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) ** 2 * sin(th0) ** 2 * sin(
                        th1) - a1 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * sin(th0) ** 3 * sin(th1) - a1 * conjugate(
                        sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * sin(th0) ** 3 * sin(th1) + a1 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(
                        th1) + a1 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(
                        th1) - a2 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(
                        th1) - a2 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) + a2 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(
                        th1) + a2 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(
                        th1) + conjugate(cos(th0)) * cos(th0) * (-(conjugate(a2 * cos(th1 + th2)) * cos(th1) * (
                        1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                            a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate(
                            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                            a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0))) + conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * cos(th1) * (1 + conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (a1 * cos(th1) + a2 * cos(
                        th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                                                               a1 * cos(th1) + a2 * cos(
                                                                                   th1 + th2)) * sin(
                        th0)) + (-(conjugate(a2 * cos(th0) * sin(th1 + th2)) * (
                        cos(th0) + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                            a1 * cos(th1) + a2 * cos(th1 + th2)))) + conjugate(
                        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                     cos(th0) + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                                         a1 * cos(th1) + a2 * cos(th1 + th2))) - (
                                     conjugate(a2 * sin(th0) * sin(th1 + th2)) - conjugate(
                                         sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2)))) * (
                                     conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                         a1 * cos(th1) + a2 * cos(th1 + th2)) + sin(th0))) * sin(th1)) - conjugate(
                        a2 * cos(th1 + th2)) * (conjugate(sin(th0)) * cos(th1) * sin(th0) * (
                        1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                            a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate(
                            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                            a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + a2 * (
                                                    conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                        cos(th0) + conjugate(
                                                            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                                                            a1 * cos(th1) + a2 * cos(th1 + th2))) + conjugate(
                                                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                        conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                                            a1 * cos(th1) + a2 * cos(th1 + th2)) + sin(th0))) * sin(
                        th2)) + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * cos(th1) * sin(th1 + th2) + a1 * a2 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) ** 2 * sin(
                        th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) ** 2 * sin(
                        th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
                        th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
                        th1 + th2) + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * sin(th0) * sin(th1 + th2) + a1 * a2 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th1) ** 2 * sin(th0) ** 2 * sin(
                        th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1) ** 2 * sin(th0) ** 2 * sin(
                        th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(
                        th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(th1 + th2)), (
                    a2 * conjugate(cos(th0)) * (
                        conjugate(a2 * cos(th1 + th2)) - conjugate(a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                        conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th0) - conjugate(
                            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * sin(th0)) * sin(th2)) / (
                    conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * cos(th1) * sin(
                        th0) + a1 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * cos(th0) * cos(
                        th1) ** 2 * sin(th0) + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * cos(th0) * cos(
                        th1) * cos(th1 + th2) * sin(th0) + a1 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) ** 2 * sin(
                        th0) ** 2 + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) * cos(th1 + th2) * sin(
                        th0) ** 2 - a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * cos(th1 + th2) * sin(th1) - a1 * a2 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
                        th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
                        th1) - a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(
                        th1) - a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(th1) - conjugate(
                        sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * sin(th0) * sin(
                        th1) + conjugate(sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(
                        th0) * sin(th0) * sin(th1) - a1 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
                        th1) - a1 * conjugate(
                        sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
                        th1) + a1 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
                        th1) + a1 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
                        th1) - a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) * sin(th1) - a2 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
                        th1) - a2 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
                        th1) + a2 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(
                        th1 + th2) * sin(th0) * sin(th1) + a2 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1 + th2) * sin(
                        th0) * sin(th1) - conjugate(sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * sin(
                        th0) ** 2 * sin(th1) + conjugate(sin(th0)) * conjugate(
                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) ** 2 * sin(
                        th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(
                        th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(
                        th1) - a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th1 + th2) ** 2 * sin(th0) ** 2 * sin(
                        th1) - a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) ** 2 * sin(th0) ** 2 * sin(
                        th1) - a1 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * sin(th0) ** 3 * sin(th1) - a1 * conjugate(
                        sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * sin(th0) ** 3 * sin(th1) + a1 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(
                        th1) + a1 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(
                        th1) - a2 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(
                        th1) - a2 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) + a2 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(
                        th1) + a2 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(
                        th1) + conjugate(cos(th0)) * cos(th0) * (-(conjugate(a2 * cos(th1 + th2)) * cos(th1) * (
                        1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                            a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate(
                            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                            a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0))) + conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * cos(th1) * (1 + conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (a1 * cos(th1) + a2 * cos(
                        th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                                                               a1 * cos(th1) + a2 * cos(
                                                                                   th1 + th2)) * sin(
                        th0)) + (-(conjugate(a2 * cos(th0) * sin(th1 + th2)) * (
                        cos(th0) + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                            a1 * cos(th1) + a2 * cos(th1 + th2)))) + conjugate(
                        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                     cos(th0) + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                                         a1 * cos(th1) + a2 * cos(th1 + th2))) - (
                                     conjugate(a2 * sin(th0) * sin(th1 + th2)) - conjugate(
                                         sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2)))) * (
                                     conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                         a1 * cos(th1) + a2 * cos(th1 + th2)) + sin(th0))) * sin(th1)) - conjugate(
                        a2 * cos(th1 + th2)) * (conjugate(sin(th0)) * cos(th1) * sin(th0) * (
                        1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                            a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate(
                            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                            a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + a2 * (
                                                    conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                        cos(th0) + conjugate(
                                                            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                                                            a1 * cos(th1) + a2 * cos(th1 + th2))) + conjugate(
                                                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                        conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                                            a1 * cos(th1) + a2 * cos(th1 + th2)) + sin(th0))) * sin(
                        th2)) + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * cos(th1) * sin(th1 + th2) + a1 * a2 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) ** 2 * sin(
                        th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) ** 2 * sin(
                        th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
                        th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
                        th1 + th2) + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * sin(th0) * sin(th1 + th2) + a1 * a2 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th1) ** 2 * sin(th0) ** 2 * sin(
                        th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1) ** 2 * sin(th0) ** 2 * sin(
                        th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(
                        th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(th1 + th2)),
                -((
                      conjugate(
                          cos(
                              th0)) * cos(
                          th0) * (
                          conjugate(
                              a2 * cos(
                                  th1 + th2)) * cos(
                              th1) - conjugate(
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2)) * cos(
                              th1) + (
                              conjugate(
                                  a2 * cos(
                                      th0) * sin(
                                      th1 + th2)) * cos(
                                  th0) - conjugate(
                                  cos(
                                      th0) * (
                                      a1 * sin(
                                          th1) + a2 * sin(
                                          th1 + th2))) * cos(
                                  th0) + (
                                  conjugate(
                                      a2 * sin(
                                          th0) * sin(
                                          th1 + th2)) - conjugate(
                                      sin(
                                          th0) * (
                                          a1 * sin(
                                              th1) + a2 * sin(
                                              th1 + th2)))) * sin(
                                  th0)) * sin(
                              th1)) + (
                          2 * conjugate(
                              sin(
                                  th0)) * sin(
                              th0) * (
                              conjugate(
                                  a2 * cos(
                                      th0) * sin(
                                      th1 + th2)) * cos(
                                  th0) - conjugate(
                                  cos(
                                      th0) * (
                                      a1 * sin(
                                          th1) + a2 * sin(
                                          th1 + th2))) * cos(
                                  th0) + (
                                  conjugate(
                                      a2 * sin(
                                          th0) * sin(
                                          th1 + th2)) - conjugate(
                                      sin(
                                          th0) * (
                                          a1 * sin(
                                              th1) + a2 * sin(
                                              th1 + th2)))) * sin(
                                  th0)) * sin(
                              th1) - 2 * conjugate(
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2)) * (
                              conjugate(
                                  sin(
                                      th0)) * cos(
                                  th1) * sin(
                                  th0) + a2 * (
                                  conjugate(
                                      a2 * cos(
                                          th0) * sin(
                                          th1 + th2)) * cos(
                                      th0) + conjugate(
                                      a2 * sin(
                                          th0) * sin(
                                          th1 + th2)) * sin(
                                      th0)) * sin(
                                  th2)) + 2 * conjugate(
                              a2 * cos(
                                  th1 + th2)) * (
                              conjugate(
                                  sin(
                                      th0)) * cos(
                                  th1) * sin(
                                  th0) + a2 * (
                                  conjugate(
                                      cos(
                                          th0) * (
                                          a1 * sin(
                                              th1) + a2 * sin(
                                              th1 + th2))) * cos(
                                      th0) + conjugate(
                                      sin(
                                          th0) * (
                                          a1 * sin(
                                              th1) + a2 * sin(
                                              th1 + th2))) * sin(
                                      th0)) * sin(
                                  th2))) / 2) / (
                      conjugate(
                          a1 * cos(
                              th1) + a2 * cos(
                              th1 + th2)) * conjugate(
                          sin(
                              th0)) * cos(
                          th1) * sin(
                          th0) + a1 * conjugate(
                          a1 * cos(
                              th1) + a2 * cos(
                              th1 + th2)) * conjugate(
                          cos(
                              th0) * (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2))) * conjugate(
                          sin(
                              th0)) * cos(
                          th0) * cos(
                          th1) ** 2 * sin(
                          th0) + a2 * conjugate(
                          a1 * cos(
                              th1) + a2 * cos(
                              th1 + th2)) * conjugate(
                          cos(
                              th0) * (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2))) * conjugate(
                          sin(
                              th0)) * cos(
                          th0) * cos(
                          th1) * cos(
                          th1 + th2) * sin(
                          th0) + a1 * conjugate(
                          a1 * cos(
                              th1) + a2 * cos(
                              th1 + th2)) * conjugate(
                          sin(
                              th0)) * conjugate(
                          (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2)) * sin(
                              th0)) * cos(
                          th1) ** 2 * sin(
                          th0) ** 2 + a2 * conjugate(
                          a1 * cos(
                              th1) + a2 * cos(
                              th1 + th2)) * conjugate(
                          sin(
                              th0)) * conjugate(
                          (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2)) * sin(
                              th0)) * cos(
                          th1) * cos(
                          th1 + th2) * sin(
                          th0) ** 2 - a2 * conjugate(
                          a1 * cos(
                              th1) + a2 * cos(
                              th1 + th2)) * conjugate(
                          a2 * cos(
                              th0) * sin(
                              th1 + th2)) * cos(
                          th0) * cos(
                          th1 + th2) * sin(
                          th1) - a1 * a2 * conjugate(
                          a1 * cos(
                              th1) + a2 * cos(
                              th1 + th2)) * conjugate(
                          cos(
                              th0) * (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2))) * conjugate(
                          a2 * cos(
                              th0) * sin(
                              th1 + th2)) * cos(
                          th0) ** 2 * cos(
                          th1) * cos(
                          th1 + th2) * sin(
                          th1) - a1 * a2 * conjugate(
                          a1 * cos(
                              th1) + a2 * cos(
                              th1 + th2)) * conjugate(
                          (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2)) * sin(
                              th0)) * conjugate(
                          a2 * sin(
                              th0) * sin(
                              th1 + th2)) * cos(
                          th0) ** 2 * cos(
                          th1) * cos(
                          th1 + th2) * sin(
                          th1) - a2 ** 2 * conjugate(
                          a1 * cos(
                              th1) + a2 * cos(
                              th1 + th2)) * conjugate(
                          cos(
                              th0) * (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2))) * conjugate(
                          a2 * cos(
                              th0) * sin(
                              th1 + th2)) * cos(
                          th0) ** 2 * cos(
                          th1 + th2) ** 2 * sin(
                          th1) - a2 ** 2 * conjugate(
                          a1 * cos(
                              th1) + a2 * cos(
                              th1 + th2)) * conjugate(
                          (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2)) * sin(
                              th0)) * conjugate(
                          a2 * sin(
                              th0) * sin(
                              th1 + th2)) * cos(
                          th0) ** 2 * cos(
                          th1 + th2) ** 2 * sin(
                          th1) - conjugate(
                          sin(
                              th0)) * conjugate(
                          a2 * cos(
                              th0) * sin(
                              th1 + th2)) * cos(
                          th0) * sin(
                          th0) * sin(
                          th1) + conjugate(
                          sin(
                              th0)) * conjugate(
                          cos(
                              th0) * (
                              a1 * sin(
                                  th1) + a2 * sin(
                                  th1 + th2))) * cos(
                          th0) * sin(
                          th0) * sin(
                          th1) - a1 * conjugate(
                          cos(
                              th0) * (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2))) * conjugate(
                          sin(
                              th0)) * conjugate(
                          a2 * cos(
                              th0) * sin(
                              th1 + th2)) * cos(
                          th0) ** 2 * cos(
                          th1) * sin(
                          th0) * sin(
                          th1) - a1 * conjugate(
                          sin(
                              th0)) * conjugate(
                          (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2)) * sin(
                              th0)) * conjugate(
                          a2 * sin(
                              th0) * sin(
                              th1 + th2)) * cos(
                          th0) ** 2 * cos(
                          th1) * sin(
                          th0) * sin(
                          th1) + a1 * conjugate(
                          cos(
                              th0) * (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2))) * conjugate(
                          sin(
                              th0)) * conjugate(
                          cos(
                              th0) * (
                              a1 * sin(
                                  th1) + a2 * sin(
                                  th1 + th2))) * cos(
                          th0) ** 2 * cos(
                          th1) * sin(
                          th0) * sin(
                          th1) + a1 * conjugate(
                          sin(
                              th0)) * conjugate(
                          (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2)) * sin(
                              th0)) * conjugate(
                          sin(
                              th0) * (
                              a1 * sin(
                                  th1) + a2 * sin(
                                  th1 + th2))) * cos(
                          th0) ** 2 * cos(
                          th1) * sin(
                          th0) * sin(
                          th1) - a2 * conjugate(
                          a1 * cos(
                              th1) + a2 * cos(
                              th1 + th2)) * conjugate(
                          a2 * sin(
                              th0) * sin(
                              th1 + th2)) * cos(
                          th1 + th2) * sin(
                          th0) * sin(
                          th1) - a2 * conjugate(
                          cos(
                              th0) * (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2))) * conjugate(
                          sin(
                              th0)) * conjugate(
                          a2 * cos(
                              th0) * sin(
                              th1 + th2)) * cos(
                          th0) ** 2 * cos(
                          th1 + th2) * sin(
                          th0) * sin(
                          th1) - a2 * conjugate(
                          sin(
                              th0)) * conjugate(
                          (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2)) * sin(
                              th0)) * conjugate(
                          a2 * sin(
                              th0) * sin(
                              th1 + th2)) * cos(
                          th0) ** 2 * cos(
                          th1 + th2) * sin(
                          th0) * sin(
                          th1) + a2 * conjugate(
                          cos(
                              th0) * (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2))) * conjugate(
                          sin(
                              th0)) * conjugate(
                          cos(
                              th0) * (
                              a1 * sin(
                                  th1) + a2 * sin(
                                  th1 + th2))) * cos(
                          th0) ** 2 * cos(
                          th1 + th2) * sin(
                          th0) * sin(
                          th1) + a2 * conjugate(
                          sin(
                              th0)) * conjugate(
                          (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2)) * sin(
                              th0)) * conjugate(
                          sin(
                              th0) * (
                              a1 * sin(
                                  th1) + a2 * sin(
                                  th1 + th2))) * cos(
                          th0) ** 2 * cos(
                          th1 + th2) * sin(
                          th0) * sin(
                          th1) - conjugate(
                          sin(
                              th0)) * conjugate(
                          a2 * sin(
                              th0) * sin(
                              th1 + th2)) * sin(
                          th0) ** 2 * sin(
                          th1) + conjugate(
                          sin(
                              th0)) * conjugate(
                          sin(
                              th0) * (
                              a1 * sin(
                                  th1) + a2 * sin(
                                  th1 + th2))) * sin(
                          th0) ** 2 * sin(
                          th1) - a1 * a2 * conjugate(
                          a1 * cos(
                              th1) + a2 * cos(
                              th1 + th2)) * conjugate(
                          cos(
                              th0) * (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2))) * conjugate(
                          a2 * cos(
                              th0) * sin(
                              th1 + th2)) * cos(
                          th1) * cos(
                          th1 + th2) * sin(
                          th0) ** 2 * sin(
                          th1) - a1 * a2 * conjugate(
                          a1 * cos(
                              th1) + a2 * cos(
                              th1 + th2)) * conjugate(
                          (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2)) * sin(
                              th0)) * conjugate(
                          a2 * sin(
                              th0) * sin(
                              th1 + th2)) * cos(
                          th1) * cos(
                          th1 + th2) * sin(
                          th0) ** 2 * sin(
                          th1) - a2 ** 2 * conjugate(
                          a1 * cos(
                              th1) + a2 * cos(
                              th1 + th2)) * conjugate(
                          cos(
                              th0) * (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2))) * conjugate(
                          a2 * cos(
                              th0) * sin(
                              th1 + th2)) * cos(
                          th1 + th2) ** 2 * sin(
                          th0) ** 2 * sin(
                          th1) - a2 ** 2 * conjugate(
                          a1 * cos(
                              th1) + a2 * cos(
                              th1 + th2)) * conjugate(
                          (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2)) * sin(
                              th0)) * conjugate(
                          a2 * sin(
                              th0) * sin(
                              th1 + th2)) * cos(
                          th1 + th2) ** 2 * sin(
                          th0) ** 2 * sin(
                          th1) - a1 * conjugate(
                          cos(
                              th0) * (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2))) * conjugate(
                          sin(
                              th0)) * conjugate(
                          a2 * cos(
                              th0) * sin(
                              th1 + th2)) * cos(
                          th1) * sin(
                          th0) ** 3 * sin(
                          th1) - a1 * conjugate(
                          sin(
                              th0)) * conjugate(
                          (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2)) * sin(
                              th0)) * conjugate(
                          a2 * sin(
                              th0) * sin(
                              th1 + th2)) * cos(
                          th1) * sin(
                          th0) ** 3 * sin(
                          th1) + a1 * conjugate(
                          cos(
                              th0) * (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2))) * conjugate(
                          sin(
                              th0)) * conjugate(
                          cos(
                              th0) * (
                              a1 * sin(
                                  th1) + a2 * sin(
                                  th1 + th2))) * cos(
                          th1) * sin(
                          th0) ** 3 * sin(
                          th1) + a1 * conjugate(
                          sin(
                              th0)) * conjugate(
                          (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2)) * sin(
                              th0)) * conjugate(
                          sin(
                              th0) * (
                              a1 * sin(
                                  th1) + a2 * sin(
                                  th1 + th2))) * cos(
                          th1) * sin(
                          th0) ** 3 * sin(
                          th1) - a2 * conjugate(
                          cos(
                              th0) * (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2))) * conjugate(
                          sin(
                              th0)) * conjugate(
                          a2 * cos(
                              th0) * sin(
                              th1 + th2)) * cos(
                          th1 + th2) * sin(
                          th0) ** 3 * sin(
                          th1) - a2 * conjugate(
                          sin(
                              th0)) * conjugate(
                          (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2)) * sin(
                              th0)) * conjugate(
                          a2 * sin(
                              th0) * sin(
                              th1 + th2)) * cos(
                          th1 + th2) * sin(
                          th0) ** 3 * sin(
                          th1) + a2 * conjugate(
                          cos(
                              th0) * (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2))) * conjugate(
                          sin(
                              th0)) * conjugate(
                          cos(
                              th0) * (
                              a1 * sin(
                                  th1) + a2 * sin(
                                  th1 + th2))) * cos(
                          th1 + th2) * sin(
                          th0) ** 3 * sin(
                          th1) + a2 * conjugate(
                          sin(
                              th0)) * conjugate(
                          (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2)) * sin(
                              th0)) * conjugate(
                          sin(
                              th0) * (
                              a1 * sin(
                                  th1) + a2 * sin(
                                  th1 + th2))) * cos(
                          th1 + th2) * sin(
                          th0) ** 3 * sin(
                          th1) + conjugate(
                          cos(
                              th0)) * cos(
                          th0) * (
                          -(
                              conjugate(
                                  a2 * cos(
                                      th1 + th2)) * cos(
                                  th1) * (
                                  1 + conjugate(
                                      cos(
                                          th0) * (
                                          a1 * cos(
                                              th1) + a2 * cos(
                                              th1 + th2))) * cos(
                                      th0) * (
                                      a1 * cos(
                                          th1) + a2 * cos(
                                          th1 + th2)) + conjugate(
                                      (
                                          a1 * cos(
                                              th1) + a2 * cos(
                                              th1 + th2)) * sin(
                                          th0)) * (
                                      a1 * cos(
                                          th1) + a2 * cos(
                                          th1 + th2)) * sin(
                                      th0))) + conjugate(
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2)) * cos(
                              th1) * (
                              1 + conjugate(
                                  cos(
                                      th0) * (
                                      a1 * cos(
                                          th1) + a2 * cos(
                                          th1 + th2))) * cos(
                                  th0) * (
                                  a1 * cos(
                                      th1) + a2 * cos(
                                      th1 + th2)) + conjugate(
                                  (
                                      a1 * cos(
                                          th1) + a2 * cos(
                                          th1 + th2)) * sin(
                                      th0)) * (
                                  a1 * cos(
                                      th1) + a2 * cos(
                                      th1 + th2)) * sin(
                                  th0)) + (
                              -(
                                  conjugate(
                                      a2 * cos(
                                          th0) * sin(
                                          th1 + th2)) * (
                                      cos(
                                          th0) + conjugate(
                                          cos(
                                              th0) * (
                                              a1 * cos(
                                                  th1) + a2 * cos(
                                                  th1 + th2))) * (
                                          a1 * cos(
                                              th1) + a2 * cos(
                                              th1 + th2)))) + conjugate(
                                  cos(
                                      th0) * (
                                      a1 * sin(
                                          th1) + a2 * sin(
                                          th1 + th2))) * (
                                  cos(
                                      th0) + conjugate(
                                      cos(
                                          th0) * (
                                          a1 * cos(
                                              th1) + a2 * cos(
                                              th1 + th2))) * (
                                      a1 * cos(
                                          th1) + a2 * cos(
                                          th1 + th2))) - (
                                  conjugate(
                                      a2 * sin(
                                          th0) * sin(
                                          th1 + th2)) - conjugate(
                                      sin(
                                          th0) * (
                                          a1 * sin(
                                              th1) + a2 * sin(
                                              th1 + th2)))) * (
                                  conjugate(
                                      (
                                          a1 * cos(
                                              th1) + a2 * cos(
                                              th1 + th2)) * sin(
                                          th0)) * (
                                      a1 * cos(
                                          th1) + a2 * cos(
                                          th1 + th2)) + sin(
                                      th0))) * sin(
                              th1)) - conjugate(
                          a2 * cos(
                              th1 + th2)) * (
                          conjugate(
                              sin(
                                  th0)) * cos(
                              th1) * sin(
                              th0) * (
                              1 + conjugate(
                                  cos(
                                      th0) * (
                                      a1 * cos(
                                          th1) + a2 * cos(
                                          th1 + th2))) * cos(
                                  th0) * (
                                  a1 * cos(
                                      th1) + a2 * cos(
                                      th1 + th2)) + conjugate(
                                  (
                                      a1 * cos(
                                          th1) + a2 * cos(
                                          th1 + th2)) * sin(
                                      th0)) * (
                                  a1 * cos(
                                      th1) + a2 * cos(
                                      th1 + th2)) * sin(
                                  th0)) + a2 * (
                              conjugate(
                                  cos(
                                      th0) * (
                                      a1 * sin(
                                          th1) + a2 * sin(
                                          th1 + th2))) * (
                                  cos(
                                      th0) + conjugate(
                                      cos(
                                          th0) * (
                                          a1 * cos(
                                              th1) + a2 * cos(
                                              th1 + th2))) * (
                                      a1 * cos(
                                          th1) + a2 * cos(
                                          th1 + th2))) + conjugate(
                                  sin(
                                      th0) * (
                                      a1 * sin(
                                          th1) + a2 * sin(
                                          th1 + th2))) * (
                                  conjugate(
                                      (
                                          a1 * cos(
                                              th1) + a2 * cos(
                                              th1 + th2)) * sin(
                                          th0)) * (
                                      a1 * cos(
                                          th1) + a2 * cos(
                                          th1 + th2)) + sin(
                                      th0))) * sin(
                              th2)) + a2 * conjugate(
                          a1 * cos(
                              th1) + a2 * cos(
                              th1 + th2)) * conjugate(
                          a2 * cos(
                              th0) * sin(
                              th1 + th2)) * cos(
                          th0) * cos(
                          th1) * sin(
                          th1 + th2) + a1 * a2 * conjugate(
                          a1 * cos(
                              th1) + a2 * cos(
                              th1 + th2)) * conjugate(
                          cos(
                              th0) * (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2))) * conjugate(
                          a2 * cos(
                              th0) * sin(
                              th1 + th2)) * cos(
                          th0) ** 2 * cos(
                          th1) ** 2 * sin(
                          th1 + th2) + a1 * a2 * conjugate(
                          a1 * cos(
                              th1) + a2 * cos(
                              th1 + th2)) * conjugate(
                          (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2)) * sin(
                              th0)) * conjugate(
                          a2 * sin(
                              th0) * sin(
                              th1 + th2)) * cos(
                          th0) ** 2 * cos(
                          th1) ** 2 * sin(
                          th1 + th2) + a2 ** 2 * conjugate(
                          a1 * cos(
                              th1) + a2 * cos(
                              th1 + th2)) * conjugate(
                          cos(
                              th0) * (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2))) * conjugate(
                          a2 * cos(
                              th0) * sin(
                              th1 + th2)) * cos(
                          th0) ** 2 * cos(
                          th1) * cos(
                          th1 + th2) * sin(
                          th1 + th2) + a2 ** 2 * conjugate(
                          a1 * cos(
                              th1) + a2 * cos(
                              th1 + th2)) * conjugate(
                          (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2)) * sin(
                              th0)) * conjugate(
                          a2 * sin(
                              th0) * sin(
                              th1 + th2)) * cos(
                          th0) ** 2 * cos(
                          th1) * cos(
                          th1 + th2) * sin(
                          th1 + th2) + a2 * conjugate(
                          a1 * cos(
                              th1) + a2 * cos(
                              th1 + th2)) * conjugate(
                          a2 * sin(
                              th0) * sin(
                              th1 + th2)) * cos(
                          th1) * sin(
                          th0) * sin(
                          th1 + th2) + a1 * a2 * conjugate(
                          a1 * cos(
                              th1) + a2 * cos(
                              th1 + th2)) * conjugate(
                          cos(
                              th0) * (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2))) * conjugate(
                          a2 * cos(
                              th0) * sin(
                              th1 + th2)) * cos(
                          th1) ** 2 * sin(
                          th0) ** 2 * sin(
                          th1 + th2) + a1 * a2 * conjugate(
                          a1 * cos(
                              th1) + a2 * cos(
                              th1 + th2)) * conjugate(
                          (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2)) * sin(
                              th0)) * conjugate(
                          a2 * sin(
                              th0) * sin(
                              th1 + th2)) * cos(
                          th1) ** 2 * sin(
                          th0) ** 2 * sin(
                          th1 + th2) + a2 ** 2 * conjugate(
                          a1 * cos(
                              th1) + a2 * cos(
                              th1 + th2)) * conjugate(
                          cos(
                              th0) * (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2))) * conjugate(
                          a2 * cos(
                              th0) * sin(
                              th1 + th2)) * cos(
                          th1) * cos(
                          th1 + th2) * sin(
                          th0) ** 2 * sin(
                          th1 + th2) + a2 ** 2 * conjugate(
                          a1 * cos(
                              th1) + a2 * cos(
                              th1 + th2)) * conjugate(
                          (
                              a1 * cos(
                                  th1) + a2 * cos(
                                  th1 + th2)) * sin(
                              th0)) * conjugate(
                          a2 * sin(
                              th0) * sin(
                              th1 + th2)) * cos(
                          th1) * cos(
                          th1 + th2) * sin(
                          th0) ** 2 * sin(
                          th1 + th2)))])

    r2 = array([-((a2 * conjugate(a2 * cos(th1 + th2)) * conjugate(
        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) + a1 * a2 * conjugate(
        a2 * cos(th1 + th2)) * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * cos(th1) * cos(th1 + th2) + a1 * a2 * conjugate(
        a2 * cos(th1 + th2)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * cos(th1) * cos(th1 + th2) + a2 ** 2 * conjugate(
        a2 * cos(th1 + th2)) * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * cos(th1 + th2) ** 2 + a2 ** 2 * conjugate(
        a2 * cos(th1 + th2)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * cos(th1 + th2) ** 2 - a2 * conjugate(
        a1 * cos(th1) + a2 * cos(th1 + th2)) * cos(th1 + th2) * (
                       conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                           a2 * sin(th0) * sin(th1 + th2)) * cos(th0) * (
                       a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate(
                           a2 * cos(th0) * sin(th1 + th2)) * (
                           1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                               a1 * cos(th1) + a2 * cos(th1 + th2)))) + conjugate(cos(th0)) * cos(th0) * (-(
        conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
            conjugate(a2 * sin(th0) * sin(th1 + th2)) - conjugate(
                sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2)))) * cos(
            th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) - conjugate(a2 * cos(th0) * sin(th1 + th2)) * (1 + conjugate(
        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (a1 * cos(th1) + a2 * cos(
        th1 + th2))) + conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (1 + conjugate(
        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (a1 * cos(th1) + a2 * cos(
        th1 + th2)))) - conjugate(sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * sin(th0) + conjugate(
        sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) - a1 * conjugate(
        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * cos(th1) * sin(th0) - a1 * conjugate(sin(th0)) * conjugate(
        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(th0) * cos(
        th1) * sin(th0) + a1 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
        sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * cos(th1) * sin(
        th0) + a1 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * cos(th1) * sin(th0) - a2 * conjugate(
        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * cos(th1 + th2) * sin(th0) - a2 * conjugate(sin(th0)) * conjugate(
        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(th0) * cos(
        th1 + th2) * sin(th0) + a2 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
        sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * cos(th1 + th2) * sin(
        th0) + a2 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * cos(th1 + th2) * sin(th0) + a2 * conjugate(
        a2 * sin(th0) * sin(th1 + th2)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) * sin(
        th1 + th2) - a2 * conjugate(a2 * cos(th0) * sin(th1 + th2)) * conjugate(
        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) * sin(th1 + th2)) / (a1 * (
        conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * cos(th1) * sin(th0) + a1 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * cos(th0) * cos(th1) ** 2 * sin(th0) + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * cos(th0) * cos(th1) * cos(th1 + th2) * sin(th0) + a1 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) ** 2 * sin(th0) ** 2 + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) * cos(th1 + th2) * sin(
            th0) ** 2 - a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * cos(
            th1 + th2) * sin(th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(th1) - a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
            th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(th1) - conjugate(
            sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * sin(th0) * sin(th1) + conjugate(
            sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * sin(th0) * sin(
            th1) - a1 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(th1) - a1 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(th1) + a1 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
            th1) + a1 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
            th1) - a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1 + th2) * sin(th0) * sin(th1) - a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) - a2 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(th1) + a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) + a2 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) - conjugate(sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * sin(th0) ** 2 * sin(
            th1) + conjugate(
            sin(th0)) * conjugate(sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) ** 2 * sin(
            th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1) - a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1 + th2) ** 2 * sin(th0) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) ** 2 * sin(th0) ** 2 * sin(th1) - a1 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * sin(th0) ** 3 * sin(th1) - a1 * conjugate(
            sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * sin(
            th0) ** 3 * sin(th1) + a1 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(
            th1) + a1 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(th1) - a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) - a2 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) + a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(
            th1) + a2 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) + conjugate(
            cos(th0)) * cos(th0) * (-(conjugate(a2 * cos(th1 + th2)) * cos(th1) * (
            1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0))) + conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * cos(th1) * (
                                        1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                                            a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate(
                                            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                            a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + (-(
            conjugate(a2 * cos(th0) * sin(th1 + th2)) * (
                cos(th0) + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                    a1 * cos(th1) + a2 * cos(th1 + th2)))) + conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                                                                    cos(
                                                                                                        th0) + conjugate(
                                                                                                        cos(th0) * (
                                                                                                            a1 * cos(
                                                                                                                th1) + a2 * cos(
                                                                                                                th1 + th2))) * (
                                                                                                        a1 * cos(
                                                                                                            th1) + a2 * cos(
                                                                                                            th1 + th2))) - (
                                                                                                conjugate(
                                                                                                    a2 * sin(th0) * sin(
                                                                                                        th1 + th2)) - conjugate(
                                                                                                    sin(th0) * (
                                                                                                    a1 * sin(
                                                                                                        th1) + a2 * sin(
                                                                                                        th1 + th2)))) * (
                                                                                                conjugate(
                                                                                                    (a1 * cos(
                                                                                                        th1) + a2 * cos(
                                                                                                        th1 + th2)) * sin(
                                                                                                        th0)) * (
                                                                                                a1 * cos(
                                                                                                    th1) + a2 * cos(
                                                                                                    th1 + th2)) + sin(
                                                                                                    th0))) * sin(
            th1)) - conjugate(a2 * cos(th1 + th2)) * (conjugate(sin(th0)) * cos(th1) * sin(th0) * (
            1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + a2 * (
                                                          conjugate(
                                                              cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                              cos(th0) + conjugate(
                                                                  cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                                                                  a1 * cos(th1) + a2 * cos(th1 + th2))) + conjugate(
                                                              sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                              conjugate(
                                                                  (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                                                  a1 * cos(th1) + a2 * cos(th1 + th2)) + sin(
                                                                  th0))) * sin(
            th2)) + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) * cos(th1) * sin(th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) ** 2 * sin(th1 + th2) + a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) ** 2 * sin(th1 + th2) + a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
            th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(th1 + th2) + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * sin(
            th0) * sin(
            th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1) ** 2 * sin(th0) ** 2 * sin(th1 + th2) + a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1) ** 2 * sin(th0) ** 2 * sin(th1 + th2) + a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(
            th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1 + th2)))), -((a2 * conjugate(a2 * cos(th1 + th2)) * conjugate(
        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) - conjugate(sin(th0)) * conjugate(
        a2 * sin(th0) * sin(th1 + th2)) * sin(th0) + conjugate(sin(th0)) * conjugate(
        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) + a1 * a2 * conjugate(
        a2 * cos(th1 + th2)) * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * cos(th1 + th2) * sin(th0) + a1 * a2 * conjugate(
        a2 * cos(th1 + th2)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * cos(th1 + th2) * sin(th0) + a2 ** 2 * conjugate(
        a2 * cos(th1 + th2)) * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) ** 2 * sin(th0) + a2 ** 2 * conjugate(
        a2 * cos(th1 + th2)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) ** 2 * sin(th0) - a1 * conjugate(
        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
        a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * sin(th0) ** 2 - a1 * conjugate(sin(th0)) * conjugate(
        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * sin(
        th0) ** 2 + a1 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 2 + a1 * conjugate(
        sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 2 - a2 * conjugate(
        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
        a2 * cos(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 2 - a2 * conjugate(sin(th0)) * conjugate(
        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
        th1 + th2) * sin(th0) ** 2 + a2 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
        sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(
        th0) ** 2 + a2 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 2 - a2 * conjugate(
        a1 * cos(th1) + a2 * cos(th1 + th2)) * cos(th1 + th2) * (conjugate(
        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * (
                                                                     a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(
        th0) + conjugate(a2 * sin(th0) * sin(th1 + th2)) * (1 + conjugate(
        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(
        th0))) + conjugate(cos(th0)) * cos(th0) * (-(conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
        conjugate(a2 * cos(th0) * sin(th1 + th2)) - conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2)))) * (
                                                         a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) - conjugate(
        a2 * sin(th0) * sin(th1 + th2)) * (1 + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
        a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + conjugate(
        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                       1 + conjugate(
                                                           (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                                           a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(
                                                           th0))) - a2 * conjugate(
        a2 * sin(th0) * sin(th1 + th2)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * sin(
        th1 + th2) + a2 * conjugate(a2 * cos(th0) * sin(th1 + th2)) * conjugate(
        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * sin(th1 + th2)) / (a1 * (
        conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * cos(th1) * sin(th0) + a1 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * cos(th0) * cos(th1) ** 2 * sin(th0) + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * cos(th0) * cos(th1) * cos(th1 + th2) * sin(th0) + a1 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) ** 2 * sin(th0) ** 2 + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) * cos(th1 + th2) * sin(
            th0) ** 2 - a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * cos(
            th1 + th2) * sin(th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(th1) - a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
            th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(th1) - conjugate(
            sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * sin(th0) * sin(th1) + conjugate(
            sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * sin(th0) * sin(
            th1) - a1 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(th1) - a1 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(th1) + a1 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
            th1) + a1 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
            th1) - a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1 + th2) * sin(th0) * sin(th1) - a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) - a2 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(th1) + a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) + a2 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) - conjugate(sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * sin(th0) ** 2 * sin(
            th1) + conjugate(
            sin(th0)) * conjugate(sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) ** 2 * sin(
            th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1) - a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1 + th2) ** 2 * sin(th0) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) ** 2 * sin(th0) ** 2 * sin(th1) - a1 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * sin(th0) ** 3 * sin(th1) - a1 * conjugate(
            sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * sin(
            th0) ** 3 * sin(th1) + a1 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(
            th1) + a1 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(th1) - a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) - a2 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) + a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(
            th1) + a2 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) + conjugate(
            cos(th0)) * cos(th0) * (-(conjugate(a2 * cos(th1 + th2)) * cos(th1) * (
            1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0))) + conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * cos(th1) * (
                                        1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                                            a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate(
                                            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                            a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + (-(
            conjugate(a2 * cos(th0) * sin(th1 + th2)) * (
                cos(th0) + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                    a1 * cos(th1) + a2 * cos(th1 + th2)))) + conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                                                                    cos(
                                                                                                        th0) + conjugate(
                                                                                                        cos(th0) * (
                                                                                                            a1 * cos(
                                                                                                                th1) + a2 * cos(
                                                                                                                th1 + th2))) * (
                                                                                                        a1 * cos(
                                                                                                            th1) + a2 * cos(
                                                                                                            th1 + th2))) - (
                                                                                                conjugate(
                                                                                                    a2 * sin(th0) * sin(
                                                                                                        th1 + th2)) - conjugate(
                                                                                                    sin(th0) * (
                                                                                                    a1 * sin(
                                                                                                        th1) + a2 * sin(
                                                                                                        th1 + th2)))) * (
                                                                                                conjugate(
                                                                                                    (a1 * cos(
                                                                                                        th1) + a2 * cos(
                                                                                                        th1 + th2)) * sin(
                                                                                                        th0)) * (
                                                                                                a1 * cos(
                                                                                                    th1) + a2 * cos(
                                                                                                    th1 + th2)) + sin(
                                                                                                    th0))) * sin(
            th1)) - conjugate(a2 * cos(th1 + th2)) * (conjugate(sin(th0)) * cos(th1) * sin(th0) * (
            1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + a2 * (
                                                          conjugate(
                                                              cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                              cos(th0) + conjugate(
                                                                  cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                                                                  a1 * cos(th1) + a2 * cos(th1 + th2))) + conjugate(
                                                              sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                              conjugate(
                                                                  (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                                                  a1 * cos(th1) + a2 * cos(th1 + th2)) + sin(
                                                                  th0))) * sin(
            th2)) + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) * cos(th1) * sin(th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) ** 2 * sin(th1 + th2) + a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) ** 2 * sin(th1 + th2) + a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
            th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(th1 + th2) + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * sin(
            th0) * sin(
            th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1) ** 2 * sin(th0) ** 2 * sin(th1 + th2) + a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1) ** 2 * sin(th0) ** 2 * sin(th1 + th2) + a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(
            th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1 + th2)))), (conjugate(cos(th0)) * (
        conjugate(a2 * cos(th1 + th2)) - conjugate(a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (1 + conjugate(
        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (a1 * cos(th1) + a2 * cos(
        th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (a1 * cos(th1) + a2 * cos(
        th1 + th2)) * sin(th0)) - conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * (conjugate(sin(th0)) * sin(th0) * (
        1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
            a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
            a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + a2 * (conjugate(a2 * cos(th0) * sin(th1 + th2)) * (
        cos(th0) + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
            a1 * cos(th1) + a2 * cos(th1 + th2))) + conjugate(a2 * sin(th0) * sin(th1 + th2)) * (conjugate(
        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (a1 * cos(th1) + a2 * cos(th1 + th2)) + sin(th0))) * sin(
        th1 + th2)) + conjugate(a2 * cos(th1 + th2)) * (conjugate(sin(th0)) * sin(th0) * (
        1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
            a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
            a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + a2 * (
                                                            conjugate(
                                                                cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                                cos(th0) + conjugate(
                                                                    cos(th0) * (
                                                                    a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                                                                    a1 * cos(th1) + a2 * cos(th1 + th2))) + conjugate(
                                                                sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                                conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(
                                                                    th0)) * (
                                                                    a1 * cos(th1) + a2 * cos(th1 + th2)) + sin(
                                                                    th0))) * sin(
        th1 + th2))) / (a1 * (
        conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * cos(th1) * sin(th0) + a1 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * cos(th0) * cos(th1) ** 2 * sin(th0) + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * cos(th0) * cos(th1) * cos(th1 + th2) * sin(th0) + a1 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) ** 2 * sin(th0) ** 2 + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) * cos(th1 + th2) * sin(
            th0) ** 2 - a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * cos(
            th1 + th2) * sin(th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(th1) - a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
            th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(th1) - conjugate(
            sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * sin(th0) * sin(th1) + conjugate(
            sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * sin(th0) * sin(
            th1) - a1 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(th1) - a1 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(th1) + a1 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
            th1) + a1 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
            th1) - a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1 + th2) * sin(th0) * sin(th1) - a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) - a2 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(th1) + a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) + a2 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) - conjugate(sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * sin(th0) ** 2 * sin(
            th1) + conjugate(
            sin(th0)) * conjugate(sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) ** 2 * sin(
            th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1) - a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1 + th2) ** 2 * sin(th0) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) ** 2 * sin(th0) ** 2 * sin(th1) - a1 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * sin(th0) ** 3 * sin(th1) - a1 * conjugate(
            sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * sin(
            th0) ** 3 * sin(th1) + a1 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(
            th1) + a1 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(th1) - a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) - a2 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) + a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(
            th1) + a2 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) + conjugate(
            cos(th0)) * cos(th0) * (-(conjugate(a2 * cos(th1 + th2)) * cos(th1) * (
            1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0))) + conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * cos(th1) * (
                                        1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                                            a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate(
                                            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                            a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + (-(
            conjugate(a2 * cos(th0) * sin(th1 + th2)) * (
                cos(th0) + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                    a1 * cos(th1) + a2 * cos(th1 + th2)))) + conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                                                                    cos(
                                                                                                        th0) + conjugate(
                                                                                                        cos(th0) * (
                                                                                                            a1 * cos(
                                                                                                                th1) + a2 * cos(
                                                                                                                th1 + th2))) * (
                                                                                                        a1 * cos(
                                                                                                            th1) + a2 * cos(
                                                                                                            th1 + th2))) - (
                                                                                                conjugate(
                                                                                                    a2 * sin(th0) * sin(
                                                                                                        th1 + th2)) - conjugate(
                                                                                                    sin(th0) * (
                                                                                                    a1 * sin(
                                                                                                        th1) + a2 * sin(
                                                                                                        th1 + th2)))) * (
                                                                                                conjugate(
                                                                                                    (a1 * cos(
                                                                                                        th1) + a2 * cos(
                                                                                                        th1 + th2)) * sin(
                                                                                                        th0)) * (
                                                                                                a1 * cos(
                                                                                                    th1) + a2 * cos(
                                                                                                    th1 + th2)) + sin(
                                                                                                    th0))) * sin(
            th1)) - conjugate(a2 * cos(th1 + th2)) * (conjugate(sin(th0)) * cos(th1) * sin(th0) * (
            1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + a2 * (
                                                          conjugate(
                                                              cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                              cos(th0) + conjugate(
                                                                  cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                                                                  a1 * cos(th1) + a2 * cos(th1 + th2))) + conjugate(
                                                              sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                              conjugate(
                                                                  (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                                                  a1 * cos(th1) + a2 * cos(th1 + th2)) + sin(
                                                                  th0))) * sin(
            th2)) + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) * cos(th1) * sin(th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) ** 2 * sin(th1 + th2) + a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) ** 2 * sin(th1 + th2) + a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
            th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(th1 + th2) + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * sin(
            th0) * sin(
            th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1) ** 2 * sin(th0) ** 2 * sin(th1 + th2) + a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1) ** 2 * sin(th0) ** 2 * sin(th1 + th2) + a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(
            th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1 + th2))), (a2 * conjugate(sin(th0)) * (-(
        conjugate(a2 * cos(th1 + th2)) * cos(th1 + th2) * (
            1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0))) + conjugate(
        a1 * cos(th1) + a2 * cos(th1 + th2)) * cos(
        th1 + th2) * (1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
        a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                          a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + (-(
        conjugate(a2 * cos(th0) * sin(th1 + th2)) * (
            cos(th0) + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)))) + conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                                                  cos(th0) + conjugate(cos(th0) * (
                                                                                      a1 * cos(th1) + a2 * cos(
                                                                                          th1 + th2))) * (
                                                                                      a1 * cos(th1) + a2 * cos(
                                                                                          th1 + th2))) - (
                                                                                  conjugate(a2 * sin(th0) * sin(
                                                                                      th1 + th2)) - conjugate(
                                                                                      sin(th0) * (
                                                                                          a1 * sin(th1) + a2 * sin(
                                                                                              th1 + th2)))) * (
                                                                                  conjugate((a1 * cos(th1) + a2 * cos(
                                                                                      th1 + th2)) * sin(th0)) * (
                                                                                      a1 * cos(th1) + a2 * cos(
                                                                                          th1 + th2)) + sin(
                                                                                      th0))) * sin(th1 + th2))) / (
                a1 * (
                    conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * cos(th1) * sin(
                        th0) + a1 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        sin(th0)) * cos(th0) * cos(th1) ** 2 * sin(th0) + a2 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        sin(th0)) * cos(th0) * cos(th1) * cos(th1 + th2) * sin(th0) + a1 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) ** 2 * sin(
                        th0) ** 2 + a2 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) * cos(th1 + th2) * sin(
                        th0) ** 2 - a2 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
                        th0) * cos(
                        th1 + th2) * sin(th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(
                        th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(th1) - a1 * a2 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
                        th1) - a2 ** 2 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(
                        th1) - a2 ** 2 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(th1) - conjugate(
                        sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * sin(th0) * sin(
                        th1) + conjugate(
                        sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * sin(
                        th0) * sin(
                        th1) - a1 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        sin(th0)) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
                        th1) - a1 * conjugate(
                        sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
                        th1) + a1 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
                        th1) + a1 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
                        th1) - a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(
                        th1 + th2) * sin(th0) * sin(th1) - a2 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(
                        th0) * sin(
                        th1) - a2 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
                        th1) + a2 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1 + th2) * sin(
                        th0) * sin(
                        th1) + a2 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1 + th2) * sin(
                        th0) * sin(
                        th1) - conjugate(sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * sin(th0) ** 2 * sin(
                        th1) + conjugate(
                        sin(th0)) * conjugate(sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) ** 2 * sin(
                        th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * cos(
                        th1 + th2) * sin(th0) ** 2 * sin(th1) - a1 * a2 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * cos(
                        th1 + th2) * sin(th0) ** 2 * sin(th1) - a2 ** 2 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(
                        th1 + th2) ** 2 * sin(th0) ** 2 * sin(th1) - a2 ** 2 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) ** 2 * sin(th0) ** 2 * sin(
                        th1) - a1 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * sin(th0) ** 3 * sin(th1) - a1 * conjugate(
                        sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * sin(
                        th0) ** 3 * sin(th1) + a1 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(
                        th0) ** 3 * sin(
                        th1) + a1 * conjugate(sin(th0)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(
                        th1) - a2 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) - a2 * conjugate(
                        sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) + a2 * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
                        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(
                        th1) + a2 * conjugate(
                        sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(
                        th1) + conjugate(
                        cos(th0)) * cos(th0) * (-(conjugate(a2 * cos(th1 + th2)) * cos(th1) * (
                        1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                            a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate(
                            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                            a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0))) + conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * cos(th1) * (
                                                    1 + conjugate(
                                                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                                                        a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate(
                                                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                                        a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + (-(
                        conjugate(a2 * cos(th0) * sin(th1 + th2)) * (
                            cos(th0) + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                                a1 * cos(th1) + a2 * cos(th1 + th2)))) + conjugate(
                        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                                                                                cos(
                                                                                                                    th0) + conjugate(
                                                                                                                    cos(
                                                                                                                        th0) * (
                                                                                                                        a1 * cos(
                                                                                                                            th1) + a2 * cos(
                                                                                                                            th1 + th2))) * (
                                                                                                                    a1 * cos(
                                                                                                                        th1) + a2 * cos(
                                                                                                                        th1 + th2))) - (
                                                                                                            conjugate(
                                                                                                                a2 * sin(
                                                                                                                    th0) * sin(
                                                                                                                    th1 + th2)) - conjugate(
                                                                                                                sin(
                                                                                                                    th0) * (
                                                                                                                a1 * sin(
                                                                                                                    th1) + a2 * sin(
                                                                                                                    th1 + th2)))) * (
                                                                                                            conjugate(
                                                                                                                (
                                                                                                                a1 * cos(
                                                                                                                    th1) + a2 * cos(
                                                                                                                    th1 + th2)) * sin(
                                                                                                                    th0)) * (
                                                                                                            a1 * cos(
                                                                                                                th1) + a2 * cos(
                                                                                                                th1 + th2)) + sin(
                                                                                                                th0))) * sin(
                        th1)) - conjugate(a2 * cos(th1 + th2)) * (conjugate(sin(th0)) * cos(th1) * sin(th0) * (
                        1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                            a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate(
                            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                            a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + a2 * (
                                                                      conjugate(cos(th0) * (
                                                                      a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                                          cos(th0) + conjugate(
                                                                              cos(th0) * (
                                                                              a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                                                                              a1 * cos(th1) + a2 * cos(
                                                                                  th1 + th2))) + conjugate(
                                                                          sin(th0) * (
                                                                          a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                                          conjugate((a1 * cos(th1) + a2 * cos(
                                                                              th1 + th2)) * sin(th0)) * (
                                                                              a1 * cos(th1) + a2 * cos(
                                                                                  th1 + th2)) + sin(th0))) * sin(
                        th2)) + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(
                        th0) * cos(th1) * sin(th1 + th2) + a1 * a2 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(
                        th0) ** 2 * cos(th1) ** 2 * sin(th1 + th2) + a1 * a2 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) ** 2 * sin(
                        th1 + th2) + a2 ** 2 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
                        th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(
                        th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(th1 + th2) + a2 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
                        th1) * sin(th0) * sin(
                        th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(
                        th1) ** 2 * sin(th0) ** 2 * sin(th1 + th2) + a1 * a2 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1) ** 2 * sin(th0) ** 2 * sin(
                        th1 + th2) + a2 ** 2 * conjugate(
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
                        a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(
                        th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
                        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                        a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * cos(
                        th1 + th2) * sin(th0) ** 2 * sin(th1 + th2))), (a2 * conjugate(cos(th0)) * (
        conjugate(a2 * cos(th1 + th2)) * cos(th1 + th2) * (
            1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) - conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * cos(
            th1 + th2) * (1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
            a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                              a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + (
        conjugate(a2 * cos(th0) * sin(th1 + th2)) * (
            cos(th0) + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                a1 * cos(th1) + a2 * cos(th1 + th2))) - conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
            cos(th0) + conjugate(cos(th0) * (
                a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                a1 * cos(th1) + a2 * cos(th1 + th2))) + (
            conjugate(a2 * sin(th0) * sin(
                th1 + th2)) - conjugate(sin(th0) * (
                a1 * sin(th1) + a2 * sin(th1 + th2)))) * (
            conjugate((a1 * cos(th1) + a2 * cos(
                th1 + th2)) * sin(th0)) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) + sin(
                th0))) * sin(th1 + th2))) / (a1 * (
        conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * cos(th1) * sin(th0) + a1 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * cos(th0) * cos(th1) ** 2 * sin(th0) + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * cos(th0) * cos(th1) * cos(th1 + th2) * sin(th0) + a1 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) ** 2 * sin(th0) ** 2 + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) * cos(th1 + th2) * sin(
            th0) ** 2 - a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * cos(
            th1 + th2) * sin(th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(th1) - a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
            th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(th1) - conjugate(
            sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * sin(th0) * sin(th1) + conjugate(
            sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * sin(th0) * sin(
            th1) - a1 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(th1) - a1 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(th1) + a1 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
            th1) + a1 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
            th1) - a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1 + th2) * sin(th0) * sin(th1) - a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) - a2 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(th1) + a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) + a2 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) - conjugate(sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * sin(th0) ** 2 * sin(
            th1) + conjugate(
            sin(th0)) * conjugate(sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) ** 2 * sin(
            th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1) - a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1 + th2) ** 2 * sin(th0) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) ** 2 * sin(th0) ** 2 * sin(th1) - a1 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * sin(th0) ** 3 * sin(th1) - a1 * conjugate(
            sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * sin(
            th0) ** 3 * sin(th1) + a1 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(
            th1) + a1 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(th1) - a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) - a2 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) + a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(
            th1) + a2 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) + conjugate(
            cos(th0)) * cos(th0) * (-(conjugate(a2 * cos(th1 + th2)) * cos(th1) * (
            1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0))) + conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * cos(th1) * (
                                        1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                                            a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate(
                                            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                            a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + (-(
            conjugate(a2 * cos(th0) * sin(th1 + th2)) * (
                cos(th0) + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                    a1 * cos(th1) + a2 * cos(th1 + th2)))) + conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                                                                    cos(
                                                                                                        th0) + conjugate(
                                                                                                        cos(th0) * (
                                                                                                            a1 * cos(
                                                                                                                th1) + a2 * cos(
                                                                                                                th1 + th2))) * (
                                                                                                        a1 * cos(
                                                                                                            th1) + a2 * cos(
                                                                                                            th1 + th2))) - (
                                                                                                conjugate(
                                                                                                    a2 * sin(th0) * sin(
                                                                                                        th1 + th2)) - conjugate(
                                                                                                    sin(th0) * (
                                                                                                    a1 * sin(
                                                                                                        th1) + a2 * sin(
                                                                                                        th1 + th2)))) * (
                                                                                                conjugate(
                                                                                                    (a1 * cos(
                                                                                                        th1) + a2 * cos(
                                                                                                        th1 + th2)) * sin(
                                                                                                        th0)) * (
                                                                                                a1 * cos(
                                                                                                    th1) + a2 * cos(
                                                                                                    th1 + th2)) + sin(
                                                                                                    th0))) * sin(
            th1)) - conjugate(a2 * cos(th1 + th2)) * (conjugate(sin(th0)) * cos(th1) * sin(th0) * (
            1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + a2 * (
                                                          conjugate(
                                                              cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                              cos(th0) + conjugate(
                                                                  cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                                                                  a1 * cos(th1) + a2 * cos(th1 + th2))) + conjugate(
                                                              sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                              conjugate(
                                                                  (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                                                  a1 * cos(th1) + a2 * cos(th1 + th2)) + sin(
                                                                  th0))) * sin(
            th2)) + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) * cos(th1) * sin(th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) ** 2 * sin(th1 + th2) + a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) ** 2 * sin(th1 + th2) + a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
            th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(th1 + th2) + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * sin(
            th0) * sin(
            th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1) ** 2 * sin(th0) ** 2 * sin(th1 + th2) + a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1) ** 2 * sin(th0) ** 2 * sin(th1 + th2) + a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(
            th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1 + th2))), -(((a1 * cos(th1) + a2 * cos(th1 + th2)) * (-(
        a2 * conjugate(a2 * cos(th1 + th2)) * conjugate(sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(
            th0) * cos(
            th1 + th2)) + conjugate(sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(th0) * sin(
        th0) - conjugate(
        sin(th0)) * conjugate(sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * sin(th0) + a2 * conjugate(
        a2 * cos(th1 + th2)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(
        th0) - conjugate(sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * sin(th0) ** 2 + conjugate(
        sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) ** 2 + a2 * conjugate(
        a1 * cos(th1) + a2 * cos(th1 + th2)) * cos(th1 + th2) * (conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
        th0) - conjugate(a2 * cos(th0) * sin(th1 + th2)) * sin(th0)) + conjugate(cos(th0)) * cos(th0) * (conjugate(
        a2 * sin(th0) * sin(th1 + th2)) * cos(th0) - conjugate(sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(
        th0) + (-conjugate(a2 * cos(th0) * sin(th1 + th2)) + conjugate(
        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2)))) * sin(th0)) + a2 * conjugate(
        a2 * sin(th0) * sin(th1 + th2)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(
        th0) ** 2 * sin(th1 + th2) - a2 * conjugate(a2 * cos(th0) * sin(th1 + th2)) * conjugate(
        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * sin(th1 + th2) + a2 * conjugate(
        a2 * sin(th0) * sin(th1 + th2)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(
        th0) ** 2 * sin(th1 + th2) - a2 * conjugate(a2 * cos(th0) * sin(th1 + th2)) * conjugate(
        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) ** 2 * sin(th1 + th2))) / (a1 * (
        conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * cos(th1) * sin(th0) + a1 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * cos(th0) * cos(th1) ** 2 * sin(th0) + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * cos(th0) * cos(th1) * cos(th1 + th2) * sin(th0) + a1 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) ** 2 * sin(th0) ** 2 + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) * cos(th1 + th2) * sin(
            th0) ** 2 - a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * cos(
            th1 + th2) * sin(th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(th1) - a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
            th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(th1) - conjugate(
            sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * sin(th0) * sin(th1) + conjugate(
            sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * sin(th0) * sin(
            th1) - a1 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(th1) - a1 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(th1) + a1 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
            th1) + a1 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
            th1) - a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1 + th2) * sin(th0) * sin(th1) - a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) - a2 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(th1) + a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) + a2 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) - conjugate(sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * sin(th0) ** 2 * sin(
            th1) + conjugate(
            sin(th0)) * conjugate(sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) ** 2 * sin(
            th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1) - a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1 + th2) ** 2 * sin(th0) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) ** 2 * sin(th0) ** 2 * sin(th1) - a1 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * sin(th0) ** 3 * sin(th1) - a1 * conjugate(
            sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * sin(
            th0) ** 3 * sin(th1) + a1 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(
            th1) + a1 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(th1) - a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) - a2 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) + a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(
            th1) + a2 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) + conjugate(
            cos(th0)) * cos(th0) * (-(conjugate(a2 * cos(th1 + th2)) * cos(th1) * (
            1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0))) + conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * cos(th1) * (
                                        1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                                            a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate(
                                            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                            a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + (-(
            conjugate(a2 * cos(th0) * sin(th1 + th2)) * (
                cos(th0) + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                    a1 * cos(th1) + a2 * cos(th1 + th2)))) + conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                                                                    cos(
                                                                                                        th0) + conjugate(
                                                                                                        cos(th0) * (
                                                                                                            a1 * cos(
                                                                                                                th1) + a2 * cos(
                                                                                                                th1 + th2))) * (
                                                                                                        a1 * cos(
                                                                                                            th1) + a2 * cos(
                                                                                                            th1 + th2))) - (
                                                                                                conjugate(
                                                                                                    a2 * sin(th0) * sin(
                                                                                                        th1 + th2)) - conjugate(
                                                                                                    sin(th0) * (
                                                                                                    a1 * sin(
                                                                                                        th1) + a2 * sin(
                                                                                                        th1 + th2)))) * (
                                                                                                conjugate(
                                                                                                    (a1 * cos(
                                                                                                        th1) + a2 * cos(
                                                                                                        th1 + th2)) * sin(
                                                                                                        th0)) * (
                                                                                                a1 * cos(
                                                                                                    th1) + a2 * cos(
                                                                                                    th1 + th2)) + sin(
                                                                                                    th0))) * sin(
            th1)) - conjugate(a2 * cos(th1 + th2)) * (conjugate(sin(th0)) * cos(th1) * sin(th0) * (
            1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + a2 * (
                                                          conjugate(
                                                              cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                              cos(th0) + conjugate(
                                                                  cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                                                                  a1 * cos(th1) + a2 * cos(th1 + th2))) + conjugate(
                                                              sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                              conjugate(
                                                                  (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                                                  a1 * cos(th1) + a2 * cos(th1 + th2)) + sin(
                                                                  th0))) * sin(
            th2)) + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) * cos(th1) * sin(th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) ** 2 * sin(th1 + th2) + a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) ** 2 * sin(th1 + th2) + a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
            th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(th1 + th2) + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * sin(
            th0) * sin(
            th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1) ** 2 * sin(th0) ** 2 * sin(th1 + th2) + a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1) ** 2 * sin(th0) ** 2 * sin(th1 + th2) + a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(
            th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1 + th2))))])

    r3 = array([-((-(
        a1 * conjugate(a2 * cos(th1 + th2)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(
            th1)) - a1 ** 2 * conjugate(a2 * cos(th1 + th2)) * conjugate(
        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * cos(th1) ** 2 - a1 ** 2 * conjugate(
        a2 * cos(th1 + th2)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * cos(th1) ** 2 - a2 * conjugate(
        a2 * cos(th1 + th2)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(
        th1 + th2) - 2 * a1 * a2 * conjugate(a2 * cos(th1 + th2)) * conjugate(
        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * cos(th1) * cos(
        th1 + th2) - 2 * a1 * a2 * conjugate(a2 * cos(th1 + th2)) * conjugate(
        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * cos(th1) * cos(th1 + th2) - a2 ** 2 * conjugate(
        a2 * cos(th1 + th2)) * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * cos(th1 + th2) ** 2 - a2 ** 2 * conjugate(
        a2 * cos(th1 + th2)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * cos(th1 + th2) ** 2 + conjugate(
        a1 * cos(th1) + a2 * cos(th1 + th2)) * (a1 * cos(th1) + a2 * cos(th1 + th2)) * (
                       conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
                           a2 * sin(th0) * sin(th1 + th2)) * cos(th0) * (
                       a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate(
                           a2 * cos(th0) * sin(th1 + th2)) * (
                           1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                               a1 * cos(th1) + a2 * cos(th1 + th2)))) + conjugate(cos(th0)) * cos(th0) * (
                       conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                           conjugate(a2 * sin(th0) * sin(th1 + th2)) - conjugate(
                               sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2)))) * cos(th0) * (
                           a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate(a2 * cos(th0) * sin(th1 + th2)) * (
                           1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                               a1 * cos(th1) + a2 * cos(th1 + th2))) - conjugate(
                           cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                           1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                               a1 * cos(th1) + a2 * cos(th1 + th2)))) + conjugate(sin(th0)) * conjugate(
        a2 * cos(th0) * sin(th1 + th2)) * sin(th0) - conjugate(sin(th0)) * conjugate(
        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) + a1 * conjugate(
        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * cos(th1) * sin(th0) + a1 * conjugate(sin(th0)) * conjugate(
        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(th0) * cos(
        th1) * sin(th0) - a1 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
        sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * cos(th1) * sin(
        th0) - a1 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * cos(th1) * sin(th0) + a2 * conjugate(
        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
        a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * cos(th1 + th2) * sin(th0) + a2 * conjugate(sin(th0)) * conjugate(
        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(th0) * cos(
        th1 + th2) * sin(th0) - a2 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
        sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * cos(th1 + th2) * sin(
        th0) - a2 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * cos(th1 + th2) * sin(th0) - a1 * conjugate(
        a2 * sin(th0) * sin(th1 + th2)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) * sin(
        th1) + a1 * conjugate(a2 * cos(th0) * sin(th1 + th2)) * conjugate(
        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) * sin(th1) - a2 * conjugate(
        a2 * sin(th0) * sin(th1 + th2)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) * sin(
        th1 + th2) + a2 * conjugate(a2 * cos(th0) * sin(th1 + th2)) * conjugate(
        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) * sin(th1 + th2)) / (a1 * (
        conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * cos(th1) * sin(th0) + a1 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * cos(th0) * cos(th1) ** 2 * sin(th0) + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * cos(th0) * cos(th1) * cos(th1 + th2) * sin(th0) + a1 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) ** 2 * sin(th0) ** 2 + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) * cos(th1 + th2) * sin(
            th0) ** 2 - a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * cos(
            th1 + th2) * sin(th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(th1) - a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
            th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(th1) - conjugate(
            sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * sin(th0) * sin(th1) + conjugate(
            sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * sin(th0) * sin(
            th1) - a1 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(th1) - a1 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(th1) + a1 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
            th1) + a1 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
            th1) - a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1 + th2) * sin(th0) * sin(th1) - a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) - a2 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(th1) + a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) + a2 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) - conjugate(sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * sin(th0) ** 2 * sin(
            th1) + conjugate(
            sin(th0)) * conjugate(sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) ** 2 * sin(
            th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1) - a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1 + th2) ** 2 * sin(th0) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) ** 2 * sin(th0) ** 2 * sin(th1) - a1 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * sin(th0) ** 3 * sin(th1) - a1 * conjugate(
            sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * sin(
            th0) ** 3 * sin(th1) + a1 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(
            th1) + a1 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(th1) - a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) - a2 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) + a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(
            th1) + a2 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) + conjugate(
            cos(th0)) * cos(th0) * (-(conjugate(a2 * cos(th1 + th2)) * cos(th1) * (
            1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0))) + conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * cos(th1) * (
                                        1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                                            a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate(
                                            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                            a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + (-(
            conjugate(a2 * cos(th0) * sin(th1 + th2)) * (
                cos(th0) + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                    a1 * cos(th1) + a2 * cos(th1 + th2)))) + conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                                                                    cos(
                                                                                                        th0) + conjugate(
                                                                                                        cos(th0) * (
                                                                                                            a1 * cos(
                                                                                                                th1) + a2 * cos(
                                                                                                                th1 + th2))) * (
                                                                                                        a1 * cos(
                                                                                                            th1) + a2 * cos(
                                                                                                            th1 + th2))) - (
                                                                                                conjugate(
                                                                                                    a2 * sin(th0) * sin(
                                                                                                        th1 + th2)) - conjugate(
                                                                                                    sin(th0) * (
                                                                                                    a1 * sin(
                                                                                                        th1) + a2 * sin(
                                                                                                        th1 + th2)))) * (
                                                                                                conjugate(
                                                                                                    (a1 * cos(
                                                                                                        th1) + a2 * cos(
                                                                                                        th1 + th2)) * sin(
                                                                                                        th0)) * (
                                                                                                a1 * cos(
                                                                                                    th1) + a2 * cos(
                                                                                                    th1 + th2)) + sin(
                                                                                                    th0))) * sin(
            th1)) - conjugate(a2 * cos(th1 + th2)) * (conjugate(sin(th0)) * cos(th1) * sin(th0) * (
            1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + a2 * (
                                                          conjugate(
                                                              cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                              cos(th0) + conjugate(
                                                                  cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                                                                  a1 * cos(th1) + a2 * cos(th1 + th2))) + conjugate(
                                                              sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                              conjugate(
                                                                  (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                                                  a1 * cos(th1) + a2 * cos(th1 + th2)) + sin(
                                                                  th0))) * sin(
            th2)) + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) * cos(th1) * sin(th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) ** 2 * sin(th1 + th2) + a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) ** 2 * sin(th1 + th2) + a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
            th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(th1 + th2) + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * sin(
            th0) * sin(
            th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1) ** 2 * sin(th0) ** 2 * sin(th1 + th2) + a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1) ** 2 * sin(th0) ** 2 * sin(th1 + th2) + a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(
            th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1 + th2)))), -((-(
        a1 * conjugate(a2 * cos(th1 + th2)) * conjugate(sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(
            th1)) - a2 * conjugate(a2 * cos(th1 + th2)) * conjugate(
        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(
        th1 + th2) + conjugate(sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * sin(th0) - conjugate(
        sin(th0)) * conjugate(sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) - a1 ** 2 * conjugate(
        a2 * cos(th1 + th2)) * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) ** 2 * sin(th0) - a1 ** 2 * conjugate(
        a2 * cos(th1 + th2)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) ** 2 * sin(th0) - 2 * a1 * a2 * conjugate(
        a2 * cos(th1 + th2)) * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * cos(th1 + th2) * sin(
        th0) - 2 * a1 * a2 * conjugate(a2 * cos(th1 + th2)) * conjugate(
        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * cos(th1 + th2) * sin(th0) - a2 ** 2 * conjugate(
        a2 * cos(th1 + th2)) * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) ** 2 * sin(th0) - a2 ** 2 * conjugate(
        a2 * cos(th1 + th2)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) ** 2 * sin(th0) + a1 * conjugate(
        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
        a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * sin(th0) ** 2 + a1 * conjugate(sin(th0)) * conjugate(
        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * sin(
        th0) ** 2 - a1 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
        cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 2 - a1 * conjugate(
        sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 2 + a2 * conjugate(
        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
        a2 * cos(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 2 + a2 * conjugate(sin(th0)) * conjugate(
        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
        th1 + th2) * sin(th0) ** 2 - a2 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
        sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(
        th0) ** 2 - a2 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 2 + conjugate(
        a1 * cos(th1) + a2 * cos(th1 + th2)) * (a1 * cos(th1) + a2 * cos(th1 + th2)) * (conjugate(
        cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * (
                                                                                            a1 * cos(th1) + a2 * cos(
                                                                                                th1 + th2)) * sin(
        th0) + conjugate(a2 * sin(th0) * sin(th1 + th2)) * (1 + conjugate(
        (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(
        th0))) + conjugate(cos(th0)) * cos(th0) * (conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
        conjugate(a2 * cos(th0) * sin(th1 + th2)) - conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2)))) * (
                                                       a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0) + conjugate(
        a2 * sin(th0) * sin(th1 + th2)) * (1 + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
        a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) - conjugate(
        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                       1 + conjugate(
                                                           (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                                           a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(
                                                           th0))) + a1 * conjugate(
        a2 * sin(th0) * sin(th1 + th2)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * sin(
        th1) - a1 * conjugate(a2 * cos(th0) * sin(th1 + th2)) * conjugate(
        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * sin(th1) + a2 * conjugate(
        a2 * sin(th0) * sin(th1 + th2)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * sin(
        th1 + th2) - a2 * conjugate(a2 * cos(th0) * sin(th1 + th2)) * conjugate(
        sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * sin(th1 + th2)) / (a1 * (
        conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * cos(th1) * sin(th0) + a1 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * cos(th0) * cos(th1) ** 2 * sin(th0) + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * cos(th0) * cos(th1) * cos(th1 + th2) * sin(th0) + a1 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) ** 2 * sin(th0) ** 2 + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) * cos(th1 + th2) * sin(
            th0) ** 2 - a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * cos(
            th1 + th2) * sin(th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(th1) - a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
            th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(th1) - conjugate(
            sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * sin(th0) * sin(th1) + conjugate(
            sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * sin(th0) * sin(
            th1) - a1 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(th1) - a1 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(th1) + a1 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
            th1) + a1 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
            th1) - a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1 + th2) * sin(th0) * sin(th1) - a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) - a2 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(th1) + a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) + a2 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) - conjugate(sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * sin(th0) ** 2 * sin(
            th1) + conjugate(
            sin(th0)) * conjugate(sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) ** 2 * sin(
            th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1) - a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1 + th2) ** 2 * sin(th0) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) ** 2 * sin(th0) ** 2 * sin(th1) - a1 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * sin(th0) ** 3 * sin(th1) - a1 * conjugate(
            sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * sin(
            th0) ** 3 * sin(th1) + a1 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(
            th1) + a1 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(th1) - a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) - a2 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) + a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(
            th1) + a2 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) + conjugate(
            cos(th0)) * cos(th0) * (-(conjugate(a2 * cos(th1 + th2)) * cos(th1) * (
            1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0))) + conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * cos(th1) * (
                                        1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                                            a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate(
                                            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                            a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + (-(
            conjugate(a2 * cos(th0) * sin(th1 + th2)) * (
                cos(th0) + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                    a1 * cos(th1) + a2 * cos(th1 + th2)))) + conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                                                                    cos(
                                                                                                        th0) + conjugate(
                                                                                                        cos(th0) * (
                                                                                                            a1 * cos(
                                                                                                                th1) + a2 * cos(
                                                                                                                th1 + th2))) * (
                                                                                                        a1 * cos(
                                                                                                            th1) + a2 * cos(
                                                                                                            th1 + th2))) - (
                                                                                                conjugate(
                                                                                                    a2 * sin(th0) * sin(
                                                                                                        th1 + th2)) - conjugate(
                                                                                                    sin(th0) * (
                                                                                                    a1 * sin(
                                                                                                        th1) + a2 * sin(
                                                                                                        th1 + th2)))) * (
                                                                                                conjugate(
                                                                                                    (a1 * cos(
                                                                                                        th1) + a2 * cos(
                                                                                                        th1 + th2)) * sin(
                                                                                                        th0)) * (
                                                                                                a1 * cos(
                                                                                                    th1) + a2 * cos(
                                                                                                    th1 + th2)) + sin(
                                                                                                    th0))) * sin(
            th1)) - conjugate(a2 * cos(th1 + th2)) * (conjugate(sin(th0)) * cos(th1) * sin(th0) * (
            1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + a2 * (
                                                          conjugate(
                                                              cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                              cos(th0) + conjugate(
                                                                  cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                                                                  a1 * cos(th1) + a2 * cos(th1 + th2))) + conjugate(
                                                              sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                              conjugate(
                                                                  (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                                                  a1 * cos(th1) + a2 * cos(th1 + th2)) + sin(
                                                                  th0))) * sin(
            th2)) + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) * cos(th1) * sin(th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) ** 2 * sin(th1 + th2) + a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) ** 2 * sin(th1 + th2) + a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
            th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(th1 + th2) + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * sin(
            th0) * sin(
            th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1) ** 2 * sin(th0) ** 2 * sin(th1 + th2) + a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1) ** 2 * sin(th0) ** 2 * sin(th1 + th2) + a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(
            th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1 + th2)))), (-(
        conjugate(cos(th0)) * (conjugate(a2 * cos(th1 + th2)) - conjugate(a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(
            th0) * (1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
            a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                        a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0))) + conjugate(
        a1 * cos(th1) + a2 * cos(th1 + th2)) * (
                                                                 conjugate(sin(th0)) * sin(th0) * (1 + conjugate(
                                                                     cos(th0) * (
                                                                     a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(
                                                                     th0) * (a1 * cos(th1) + a2 * cos(
                                                                     th1 + th2)) + conjugate(
                                                                     (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(
                                                                         th0)) * (
                                                                                                       a1 * cos(
                                                                                                           th1) + a2 * cos(
                                                                                                           th1 + th2)) * sin(
                                                                     th0)) + (
                                                                 conjugate(a2 * cos(th0) * sin(th1 + th2)) * (
                                                                     cos(th0) + conjugate(
                                                                         cos(th0) * (
                                                                         a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                                                                         a1 * cos(th1) + a2 * cos(
                                                                             th1 + th2))) + conjugate(
                                                                     a2 * sin(th0) * sin(th1 + th2)) * (conjugate(
                                                                     (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(
                                                                         th0)) * (
                                                                                                            a1 * cos(
                                                                                                                th1) + a2 * cos(
                                                                                                                th1 + th2)) + sin(
                                                                     th0))) * (
                                                                     a1 * sin(th1) + a2 * sin(th1 + th2))) - conjugate(
        a2 * cos(th1 + th2)) * (conjugate(sin(th0)) * sin(th0) * (
        1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
            a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
            a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + (
                                conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                    cos(th0) + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                                        a1 * cos(th1) + a2 * cos(th1 + th2))) + conjugate(
                                    sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                    conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                        a1 * cos(th1) + a2 * cos(th1 + th2)) + sin(th0))) * (
                                    a1 * sin(th1) + a2 * sin(th1 + th2)))) / (a1 * (
        conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * cos(th1) * sin(th0) + a1 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * cos(th0) * cos(th1) ** 2 * sin(th0) + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * cos(th0) * cos(th1) * cos(th1 + th2) * sin(th0) + a1 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) ** 2 * sin(th0) ** 2 + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) * cos(th1 + th2) * sin(
            th0) ** 2 - a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * cos(
            th1 + th2) * sin(th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(th1) - a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
            th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(th1) - conjugate(
            sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * sin(th0) * sin(th1) + conjugate(
            sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * sin(th0) * sin(
            th1) - a1 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(th1) - a1 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(th1) + a1 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
            th1) + a1 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
            th1) - a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1 + th2) * sin(th0) * sin(th1) - a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) - a2 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(th1) + a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) + a2 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) - conjugate(sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * sin(th0) ** 2 * sin(
            th1) + conjugate(
            sin(th0)) * conjugate(sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) ** 2 * sin(
            th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1) - a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1 + th2) ** 2 * sin(th0) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) ** 2 * sin(th0) ** 2 * sin(th1) - a1 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * sin(th0) ** 3 * sin(th1) - a1 * conjugate(
            sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * sin(
            th0) ** 3 * sin(th1) + a1 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(
            th1) + a1 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(th1) - a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) - a2 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) + a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(
            th1) + a2 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) + conjugate(
            cos(th0)) * cos(th0) * (-(conjugate(a2 * cos(th1 + th2)) * cos(th1) * (
            1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0))) + conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * cos(th1) * (
                                        1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                                            a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate(
                                            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                            a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + (-(
            conjugate(a2 * cos(th0) * sin(th1 + th2)) * (
                cos(th0) + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                    a1 * cos(th1) + a2 * cos(th1 + th2)))) + conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                                                                    cos(
                                                                                                        th0) + conjugate(
                                                                                                        cos(th0) * (
                                                                                                            a1 * cos(
                                                                                                                th1) + a2 * cos(
                                                                                                                th1 + th2))) * (
                                                                                                        a1 * cos(
                                                                                                            th1) + a2 * cos(
                                                                                                            th1 + th2))) - (
                                                                                                conjugate(
                                                                                                    a2 * sin(th0) * sin(
                                                                                                        th1 + th2)) - conjugate(
                                                                                                    sin(th0) * (
                                                                                                    a1 * sin(
                                                                                                        th1) + a2 * sin(
                                                                                                        th1 + th2)))) * (
                                                                                                conjugate(
                                                                                                    (a1 * cos(
                                                                                                        th1) + a2 * cos(
                                                                                                        th1 + th2)) * sin(
                                                                                                        th0)) * (
                                                                                                a1 * cos(
                                                                                                    th1) + a2 * cos(
                                                                                                    th1 + th2)) + sin(
                                                                                                    th0))) * sin(
            th1)) - conjugate(a2 * cos(th1 + th2)) * (conjugate(sin(th0)) * cos(th1) * sin(th0) * (
            1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + a2 * (
                                                          conjugate(
                                                              cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                              cos(th0) + conjugate(
                                                                  cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                                                                  a1 * cos(th1) + a2 * cos(th1 + th2))) + conjugate(
                                                              sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                              conjugate(
                                                                  (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                                                  a1 * cos(th1) + a2 * cos(th1 + th2)) + sin(
                                                                  th0))) * sin(
            th2)) + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) * cos(th1) * sin(th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) ** 2 * sin(th1 + th2) + a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) ** 2 * sin(th1 + th2) + a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
            th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(th1 + th2) + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * sin(
            th0) * sin(
            th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1) ** 2 * sin(th0) ** 2 * sin(th1 + th2) + a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1) ** 2 * sin(th0) ** 2 * sin(th1 + th2) + a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(
            th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1 + th2))), (conjugate(sin(th0)) * (
        conjugate(a2 * cos(th1 + th2)) * (a1 * cos(th1) + a2 * cos(th1 + th2)) * (
            1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) - conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * (
            a1 * cos(th1) + a2 * cos(th1 + th2)) * (
            1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + (conjugate(a2 * cos(th0) * sin(th1 + th2)) * (
            cos(th0) + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                a1 * cos(th1) + a2 * cos(th1 + th2))) - conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                                        cos(th0) + conjugate(
                                                                            cos(th0) * (
                                                                            a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                                                                            a1 * cos(th1) + a2 * cos(th1 + th2))) + (
                                                                        conjugate(
                                                                            a2 * sin(th0) * sin(th1 + th2)) - conjugate(
                                                                            sin(th0) * (
                                                                            a1 * sin(th1) + a2 * sin(th1 + th2)))) * (
                                                                        conjugate(
                                                                            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(
                                                                                th0)) * (
                                                                            a1 * cos(th1) + a2 * cos(th1 + th2)) + sin(
                                                                            th0))) * (
            a1 * sin(th1) + a2 * sin(th1 + th2)))) / (a1 * (
        conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * cos(th1) * sin(th0) + a1 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * cos(th0) * cos(th1) ** 2 * sin(th0) + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * cos(th0) * cos(th1) * cos(th1 + th2) * sin(th0) + a1 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) ** 2 * sin(th0) ** 2 + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) * cos(th1 + th2) * sin(
            th0) ** 2 - a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * cos(
            th1 + th2) * sin(th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(th1) - a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
            th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(th1) - conjugate(
            sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * sin(th0) * sin(th1) + conjugate(
            sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * sin(th0) * sin(
            th1) - a1 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(th1) - a1 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(th1) + a1 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
            th1) + a1 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
            th1) - a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1 + th2) * sin(th0) * sin(th1) - a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) - a2 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(th1) + a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) + a2 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) - conjugate(sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * sin(th0) ** 2 * sin(
            th1) + conjugate(
            sin(th0)) * conjugate(sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) ** 2 * sin(
            th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1) - a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1 + th2) ** 2 * sin(th0) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) ** 2 * sin(th0) ** 2 * sin(th1) - a1 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * sin(th0) ** 3 * sin(th1) - a1 * conjugate(
            sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * sin(
            th0) ** 3 * sin(th1) + a1 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(
            th1) + a1 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(th1) - a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) - a2 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) + a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(
            th1) + a2 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) + conjugate(
            cos(th0)) * cos(th0) * (-(conjugate(a2 * cos(th1 + th2)) * cos(th1) * (
            1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0))) + conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * cos(th1) * (
                                        1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                                            a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate(
                                            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                            a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + (-(
            conjugate(a2 * cos(th0) * sin(th1 + th2)) * (
                cos(th0) + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                    a1 * cos(th1) + a2 * cos(th1 + th2)))) + conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                                                                    cos(
                                                                                                        th0) + conjugate(
                                                                                                        cos(th0) * (
                                                                                                            a1 * cos(
                                                                                                                th1) + a2 * cos(
                                                                                                                th1 + th2))) * (
                                                                                                        a1 * cos(
                                                                                                            th1) + a2 * cos(
                                                                                                            th1 + th2))) - (
                                                                                                conjugate(
                                                                                                    a2 * sin(th0) * sin(
                                                                                                        th1 + th2)) - conjugate(
                                                                                                    sin(th0) * (
                                                                                                    a1 * sin(
                                                                                                        th1) + a2 * sin(
                                                                                                        th1 + th2)))) * (
                                                                                                conjugate(
                                                                                                    (a1 * cos(
                                                                                                        th1) + a2 * cos(
                                                                                                        th1 + th2)) * sin(
                                                                                                        th0)) * (
                                                                                                a1 * cos(
                                                                                                    th1) + a2 * cos(
                                                                                                    th1 + th2)) + sin(
                                                                                                    th0))) * sin(
            th1)) - conjugate(a2 * cos(th1 + th2)) * (conjugate(sin(th0)) * cos(th1) * sin(th0) * (
            1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + a2 * (
                                                          conjugate(
                                                              cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                              cos(th0) + conjugate(
                                                                  cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                                                                  a1 * cos(th1) + a2 * cos(th1 + th2))) + conjugate(
                                                              sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                              conjugate(
                                                                  (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                                                  a1 * cos(th1) + a2 * cos(th1 + th2)) + sin(
                                                                  th0))) * sin(
            th2)) + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) * cos(th1) * sin(th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) ** 2 * sin(th1 + th2) + a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) ** 2 * sin(th1 + th2) + a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
            th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(th1 + th2) + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * sin(
            th0) * sin(
            th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1) ** 2 * sin(th0) ** 2 * sin(th1 + th2) + a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1) ** 2 * sin(th0) ** 2 * sin(th1 + th2) + a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(
            th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1 + th2))), (conjugate(cos(th0)) * (-(
        conjugate(a2 * cos(th1 + th2)) * (a1 * cos(th1) + a2 * cos(th1 + th2)) * (
            1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0))) + conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * (
                                                                                       a1 * cos(th1) + a2 * cos(
                                                                                           th1 + th2)) * (
                                                                                       1 + conjugate(cos(th0) * (
                                                                                           a1 * cos(th1) + a2 * cos(
                                                                                               th1 + th2))) * cos(
                                                                                           th0) * (
                                                                                           a1 * cos(th1) + a2 * cos(
                                                                                               th1 + th2)) + conjugate(
                                                                                           (a1 * cos(
                                                                                               th1) + a2 * cos(
                                                                                               th1 + th2)) * sin(
                                                                                               th0)) * (
                                                                                       a1 * cos(th1) + a2 * cos(
                                                                                           th1 + th2)) * sin(th0)) - (
                                                                                   conjugate(
                                                                                       a2 * cos(th0) * sin(
                                                                                           th1 + th2)) * (
                                                                                   cos(th0) + conjugate(cos(th0) * (
                                                                                   a1 * cos(th1) + a2 * cos(
                                                                                       th1 + th2))) * (
                                                                                       a1 * cos(th1) + a2 * cos(
                                                                                           th1 + th2))) - conjugate(
                                                                                       cos(th0) * (
                                                                                       a1 * sin(th1) + a2 * sin(
                                                                                           th1 + th2))) * (cos(
                                                                                       th0) + conjugate(cos(th0) * (
                                                                                   a1 * cos(th1) + a2 * cos(
                                                                                       th1 + th2))) * (a1 * cos(
                                                                                       th1) + a2 * cos(th1 + th2))) + (
                                                                                       conjugate(
                                                                                           a2 * sin(
                                                                                               th0) * sin(
                                                                                               th1 + th2)) - conjugate(
                                                                                           sin(
                                                                                               th0) * (
                                                                                               a1 * sin(
                                                                                                   th1) + a2 * sin(
                                                                                                   th1 + th2)))) * (
                                                                                       conjugate(
                                                                                           (
                                                                                               a1 * cos(
                                                                                                   th1) + a2 * cos(
                                                                                                   th1 + th2)) * sin(
                                                                                               th0)) * (
                                                                                           a1 * cos(
                                                                                               th1) + a2 * cos(
                                                                                               th1 + th2)) + sin(
                                                                                           th0))) * (
                                                                                       a1 * sin(th1) + a2 * sin(
                                                                                           th1 + th2)))) / (a1 * (
        conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * cos(th1) * sin(th0) + a1 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * cos(th0) * cos(th1) ** 2 * sin(th0) + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * cos(th0) * cos(th1) * cos(th1 + th2) * sin(th0) + a1 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) ** 2 * sin(th0) ** 2 + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) * cos(th1 + th2) * sin(
            th0) ** 2 - a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * cos(
            th1 + th2) * sin(th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(th1) - a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
            th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(th1) - conjugate(
            sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * sin(th0) * sin(th1) + conjugate(
            sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * sin(th0) * sin(
            th1) - a1 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(th1) - a1 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(th1) + a1 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
            th1) + a1 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
            th1) - a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1 + th2) * sin(th0) * sin(th1) - a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) - a2 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(th1) + a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) + a2 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) - conjugate(sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * sin(th0) ** 2 * sin(
            th1) + conjugate(
            sin(th0)) * conjugate(sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) ** 2 * sin(
            th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1) - a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1 + th2) ** 2 * sin(th0) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) ** 2 * sin(th0) ** 2 * sin(th1) - a1 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * sin(th0) ** 3 * sin(th1) - a1 * conjugate(
            sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * sin(
            th0) ** 3 * sin(th1) + a1 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(
            th1) + a1 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(th1) - a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) - a2 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) + a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(
            th1) + a2 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) + conjugate(
            cos(th0)) * cos(th0) * (-(conjugate(a2 * cos(th1 + th2)) * cos(th1) * (
            1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0))) + conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * cos(th1) * (
                                        1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                                            a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate(
                                            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                            a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + (-(
            conjugate(a2 * cos(th0) * sin(th1 + th2)) * (
                cos(th0) + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                    a1 * cos(th1) + a2 * cos(th1 + th2)))) + conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                                                                    cos(
                                                                                                        th0) + conjugate(
                                                                                                        cos(th0) * (
                                                                                                            a1 * cos(
                                                                                                                th1) + a2 * cos(
                                                                                                                th1 + th2))) * (
                                                                                                        a1 * cos(
                                                                                                            th1) + a2 * cos(
                                                                                                            th1 + th2))) - (
                                                                                                conjugate(
                                                                                                    a2 * sin(th0) * sin(
                                                                                                        th1 + th2)) - conjugate(
                                                                                                    sin(th0) * (
                                                                                                    a1 * sin(
                                                                                                        th1) + a2 * sin(
                                                                                                        th1 + th2)))) * (
                                                                                                conjugate(
                                                                                                    (a1 * cos(
                                                                                                        th1) + a2 * cos(
                                                                                                        th1 + th2)) * sin(
                                                                                                        th0)) * (
                                                                                                a1 * cos(
                                                                                                    th1) + a2 * cos(
                                                                                                    th1 + th2)) + sin(
                                                                                                    th0))) * sin(
            th1)) - conjugate(a2 * cos(th1 + th2)) * (conjugate(sin(th0)) * cos(th1) * sin(th0) * (
            1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + a2 * (
                                                          conjugate(
                                                              cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                              cos(th0) + conjugate(
                                                                  cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                                                                  a1 * cos(th1) + a2 * cos(th1 + th2))) + conjugate(
                                                              sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                              conjugate(
                                                                  (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                                                  a1 * cos(th1) + a2 * cos(th1 + th2)) + sin(
                                                                  th0))) * sin(
            th2)) + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) * cos(th1) * sin(th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) ** 2 * sin(th1 + th2) + a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) ** 2 * sin(th1 + th2) + a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
            th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(th1 + th2) + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * sin(
            th0) * sin(
            th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1) ** 2 * sin(th0) ** 2 * sin(th1 + th2) + a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1) ** 2 * sin(th0) ** 2 * sin(th1 + th2) + a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(
            th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1 + th2))), -(((a1 * cos(th1) + a2 * cos(th1 + th2)) * (
        a1 * conjugate(a2 * cos(th1 + th2)) * conjugate(sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(
            th0) * cos(
            th1) + a2 * conjugate(a2 * cos(th1 + th2)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(
            th0) * cos(th1 + th2) - conjugate(sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(th0) * sin(
            th0) + conjugate(sin(th0)) * conjugate(sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * sin(
            th0) - a1 * conjugate(a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(
            th1) * sin(th0) - a2 * conjugate(a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) + conjugate(
            sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * sin(th0) ** 2 - conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) ** 2 - conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * (a1 * cos(th1) + a2 * cos(th1 + th2)) * (
            conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(th0) - conjugate(a2 * cos(th0) * sin(th1 + th2)) * sin(
                th0)) - conjugate(cos(th0)) * cos(th0) * (
        conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(th0) - conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) + (
            -conjugate(a2 * cos(th0) * sin(th1 + th2)) + conjugate(
                cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2)))) * sin(
            th0)) - a1 * conjugate(a2 * sin(th0) * sin(th1 + th2)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * sin(th1) + a1 * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * conjugate(sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(
            th0) ** 2 * sin(th1) - a1 * conjugate(a2 * sin(th0) * sin(th1 + th2)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) ** 2 * sin(th1) + a1 * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * conjugate(sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(
            th0) ** 2 * sin(th1) - a2 * conjugate(a2 * sin(th0) * sin(th1 + th2)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * sin(th1 + th2) + a2 * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * conjugate(sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(
            th0) ** 2 * sin(th1 + th2) - a2 * conjugate(a2 * sin(th0) * sin(th1 + th2)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) ** 2 * sin(th1 + th2) + a2 * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * conjugate(sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(
            th0) ** 2 * sin(th1 + th2))) / (a1 * (
        conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * cos(th1) * sin(th0) + a1 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * cos(th0) * cos(th1) ** 2 * sin(th0) + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * cos(th0) * cos(th1) * cos(th1 + th2) * sin(th0) + a1 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) ** 2 * sin(th0) ** 2 + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * cos(th1) * cos(th1 + th2) * sin(
            th0) ** 2 - a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * cos(
            th1 + th2) * sin(th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(th1) - a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
            th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) ** 2 * sin(th1) - conjugate(
            sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) * sin(th0) * sin(th1) + conjugate(
            sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) * sin(th0) * sin(
            th1) - a1 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(th1) - a1 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(th1) + a1 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
            th1) + a1 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1) * sin(th0) * sin(
            th1) - a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1 + th2) * sin(th0) * sin(th1) - a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) - a2 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(th1) + a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) + a2 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th0) ** 2 * cos(th1 + th2) * sin(th0) * sin(
            th1) - conjugate(sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * sin(th0) ** 2 * sin(
            th1) + conjugate(
            sin(th0)) * conjugate(sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * sin(th0) ** 2 * sin(
            th1) - a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1) - a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1 + th2) ** 2 * sin(th0) ** 2 * sin(th1) - a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) ** 2 * sin(th0) ** 2 * sin(th1) - a1 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * sin(th0) ** 3 * sin(th1) - a1 * conjugate(
            sin(th0)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * sin(
            th0) ** 3 * sin(th1) + a1 * conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            sin(th0)) * conjugate(cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(
            th1) + a1 * conjugate(sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1) * sin(th0) ** 3 * sin(th1) - a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) - a2 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) + a2 * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(sin(th0)) * conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(
            th1) + a2 * conjugate(
            sin(th0)) * conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * cos(th1 + th2) * sin(th0) ** 3 * sin(th1) + conjugate(
            cos(th0)) * cos(th0) * (-(conjugate(a2 * cos(th1 + th2)) * cos(th1) * (
            1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0))) + conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * cos(th1) * (
                                        1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                                            a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate(
                                            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                            a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + (-(
            conjugate(a2 * cos(th0) * sin(th1 + th2)) * (
                cos(th0) + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                    a1 * cos(th1) + a2 * cos(th1 + th2)))) + conjugate(
            cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                                                                    cos(
                                                                                                        th0) + conjugate(
                                                                                                        cos(th0) * (
                                                                                                            a1 * cos(
                                                                                                                th1) + a2 * cos(
                                                                                                                th1 + th2))) * (
                                                                                                        a1 * cos(
                                                                                                            th1) + a2 * cos(
                                                                                                            th1 + th2))) - (
                                                                                                conjugate(
                                                                                                    a2 * sin(th0) * sin(
                                                                                                        th1 + th2)) - conjugate(
                                                                                                    sin(th0) * (
                                                                                                    a1 * sin(
                                                                                                        th1) + a2 * sin(
                                                                                                        th1 + th2)))) * (
                                                                                                conjugate(
                                                                                                    (a1 * cos(
                                                                                                        th1) + a2 * cos(
                                                                                                        th1 + th2)) * sin(
                                                                                                        th0)) * (
                                                                                                a1 * cos(
                                                                                                    th1) + a2 * cos(
                                                                                                    th1 + th2)) + sin(
                                                                                                    th0))) * sin(
            th1)) - conjugate(a2 * cos(th1 + th2)) * (conjugate(sin(th0)) * cos(th1) * sin(th0) * (
            1 + conjugate(cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * cos(th0) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) + conjugate((a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) + a2 * (
                                                          conjugate(
                                                              cos(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                              cos(th0) + conjugate(
                                                                  cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * (
                                                                  a1 * cos(th1) + a2 * cos(th1 + th2))) + conjugate(
                                                              sin(th0) * (a1 * sin(th1) + a2 * sin(th1 + th2))) * (
                                                              conjugate(
                                                                  (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * (
                                                                  a1 * cos(th1) + a2 * cos(th1 + th2)) + sin(
                                                                  th0))) * sin(
            th2)) + a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) * cos(th1) * sin(th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) ** 2 * sin(th1 + th2) + a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) ** 2 * sin(th1 + th2) + a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(
            th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th0) ** 2 * cos(th1) * cos(th1 + th2) * sin(th1 + th2) + a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(th1) * sin(
            th0) * sin(
            th1 + th2) + a1 * a2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(a2 * cos(th0) * sin(th1 + th2)) * cos(
            th1) ** 2 * sin(th0) ** 2 * sin(th1 + th2) + a1 * a2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(
            a2 * sin(th0) * sin(th1 + th2)) * cos(th1) ** 2 * sin(th0) ** 2 * sin(th1 + th2) + a2 ** 2 * conjugate(
            a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            cos(th0) * (a1 * cos(th1) + a2 * cos(th1 + th2))) * conjugate(
            a2 * cos(th0) * sin(th1 + th2)) * cos(th1) * cos(th1 + th2) * sin(th0) ** 2 * sin(
            th1 + th2) + a2 ** 2 * conjugate(a1 * cos(th1) + a2 * cos(th1 + th2)) * conjugate(
            (a1 * cos(th1) + a2 * cos(th1 + th2)) * sin(th0)) * conjugate(a2 * sin(th0) * sin(th1 + th2)) * cos(
            th1) * cos(
            th1 + th2) * sin(th0) ** 2 * sin(th1 + th2))))])

    return array([r1, r2, r3])
