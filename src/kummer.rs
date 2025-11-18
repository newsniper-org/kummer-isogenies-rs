use crate::coding::{UnitMessageEncodable, encode};
use crate::field_types::*;
use crate::params::PARAMS;

pub fn hadamard(p: &KummerPoint) -> KummerPoint {
    let t1 = p[0].add_mod(&p[1], &PARAMS.p);
    let t2 = p[2].add_mod(&p[3], &PARAMS.p);
    let t3 = p[0].sub_mod(&p[1], &PARAMS.p);
    let t4 = p[2].sub_mod(&p[3], &PARAMS.p);
    [
        t1.add_mod(&t2, &PARAMS.p),
        t1.sub_mod(&t2, &PARAMS.p),
        t3.add_mod(&t4, &PARAMS.p),
        t3.sub_mod(&t4, &PARAMS.p),
    ]
}

pub fn squaring(p: &KummerPoint) -> KummerPoint {
    p.map(|f| f.square_mod(&PARAMS.p, PARAMS.mu))
}

pub fn four_way_mult(lhs: &KummerPoint, rhs: &KummerPoint) -> KummerPoint {
    [0usize, 1, 2, 3].map(|i| lhs[i].mul_mod(&rhs[i], &PARAMS.p, PARAMS.mu))
}

pub fn invert_4_constants(p: &KummerPoint) -> KummerPoint {
    let (mut pi1, mut pi3) = (
        p[2].mul_mod(&p[3], &PARAMS.p, PARAMS.mu),
        p[0].mul_mod(&p[1], &PARAMS.p, PARAMS.mu)
    );
    let (pi2, pi4) = (
        pi1.mul_mod(&p[0], &PARAMS.p, PARAMS.mu),
        pi3.mul_mod(&p[2], &PARAMS.p, PARAMS.mu),
    );
    (pi1, pi3) = (
        pi1.mul_mod(&p[1], &PARAMS.p, PARAMS.mu),
        pi3.mul_mod(&p[3], &PARAMS.p, PARAMS.mu)
    );

    [pi1, pi2, pi3, pi4]
}

pub fn is_on_kummer(p: &KummerPoint, fund_thetas: &[Fp2; 4]) -> bool {
    let [a, b, c, d] = *fund_thetas;
    let fund_theta_squares = squaring(fund_thetas);
    let fund_theta_quads = squaring(&fund_theta_squares);
    let [a2, b2, c2, d2] = fund_theta_squares;
    let [a4, b4, c4, d4] = fund_theta_quads;

    let [aa2, bb2, cc2, dd2] = hadamard(&fund_theta_squares);
    let p1 = a.mul_mod(&b, &PARAMS.p, PARAMS.mu).mul_mod(&c, &PARAMS.p, PARAMS.mu).mul_mod(&d, &PARAMS.p, PARAMS.mu);
    let p2 = aa2.mul_mod(&bb2, &PARAMS.p, PARAMS.mu).mul_mod(&cc2, &PARAMS.p, PARAMS.mu).mul_mod(&dd2, &PARAMS.p, PARAMS.mu);

    let mm = p1.mul_mod(&p2, &PARAMS.p, PARAMS.mu);
    
    let a2b2 = a2.mul_mod(&b2, &PARAMS.p, PARAMS.mu);
    let a2c2 = a2.mul_mod(&c2, &PARAMS.p, PARAMS.mu);
    let a2d2 = a2.mul_mod(&d2, &PARAMS.p, PARAMS.mu);
    let b2c2 = b2.mul_mod(&c2, &PARAMS.p, PARAMS.mu);
    let b2d2 = b2.mul_mod(&d2, &PARAMS.p, PARAMS.mu);
    let c2d2 = c2.mul_mod(&d2, &PARAMS.p, PARAMS.mu);

    let (m1, m2, m3) = (
        a2d2.sub_mod(&b2c2, &PARAMS.p).inv_mod(&PARAMS.p, PARAMS.mu).unwrap(),
        b2d2.sub_mod(&a2c2, &PARAMS.p).inv_mod(&PARAMS.p, PARAMS.mu).unwrap(),
        c2d2.sub_mod(&a2b2, &PARAMS.p).inv_mod(&PARAMS.p, PARAMS.mu).unwrap()
    );
    
    let ee = mm.mul_mod(&m1, &PARAMS.p, PARAMS.mu).mul_mod(&m2, &PARAMS.p, PARAMS.mu).mul_mod(&m3, &PARAMS.p, PARAMS.mu);
    let ff = m1.mul_mod(&a4.add_mod(&d4, &PARAMS.p).sub_mod(&b4, &PARAMS.p).sub_mod(&c4, &PARAMS.p), &PARAMS.p, PARAMS.mu);
    let gg = m2.mul_mod(&b4.add_mod(&d4, &PARAMS.p).sub_mod(&a4, &PARAMS.p).sub_mod(&c4, &PARAMS.p), &PARAMS.p, PARAMS.mu);
    let hh = m3.mul_mod(&c4.add_mod(&d4, &PARAMS.p).sub_mod(&a4, &PARAMS.p).sub_mod(&b4, &PARAMS.p), &PARAMS.p, PARAMS.mu);

    let [x, y, z, t] = *p;
    let ([x2, y2, z2, t2], [x4, y4, z4, t4]) = {
        let interm = squaring(p);
        (interm, squaring(&interm))
    };

    let two = Fp2::new(Fp::from_limb(2u64), Fp::zero());

    let x2y2 = x2.mul_mod(&y2, &PARAMS.p, PARAMS.mu);
    let x2z2 = x2.mul_mod(&z2, &PARAMS.p, PARAMS.mu);
    let x2t2 = x2.mul_mod(&t2, &PARAMS.p, PARAMS.mu);
    let y2z2 = y2.mul_mod(&z2, &PARAMS.p, PARAMS.mu);
    let y2t2 = y2.mul_mod(&t2, &PARAMS.p, PARAMS.mu);
    let z2t2 = z2.mul_mod(&t2, &PARAMS.p, PARAMS.mu);
    


    let tt1 = x4.add_mod(&y4, &PARAMS.p).add_mod(&z4, &PARAMS.p).add_mod(&t4, &PARAMS.p);
    let tt2 = ee.mul_mod(&x, &PARAMS.p, PARAMS.mu).mul_mod(&y, &PARAMS.p, PARAMS.mu).mul_mod(&z, &PARAMS.p, PARAMS.mu).mul_mod(&t, &PARAMS.p, PARAMS.mu).mul_mod(&two, &PARAMS.p, PARAMS.mu);
    let (tt3, tt4, tt5) = (
        x2t2.add_mod(&y2z2, &PARAMS.p).mul_mod(&ff, &PARAMS.p, PARAMS.mu),
        x2z2.add_mod(&y2t2, &PARAMS.p).mul_mod(&gg, &PARAMS.p, PARAMS.mu),
        x2y2.add_mod(&z2t2, &PARAMS.p).mul_mod(&hh, &PARAMS.p, PARAMS.mu),
    );

    tt1.add_mod(&tt2, &PARAMS.p) == tt3.add_mod(&tt4, &PARAMS.p).add_mod(&tt5, &PARAMS.p)
}

pub fn get_triple_constants_from_fund_thetas(fund_thetas: &[Fp2; 4]) -> [[Fp2; 4]; 5] {
    let s = squaring(&fund_thetas);
    let hs = hadamard(&s);
    let k3 = invert_4_constants(fund_thetas);
    let k4 = invert_4_constants(&hs);
    [*fund_thetas, s, hs, k3, k4]
}

pub fn normalise(p: &KummerPoint) -> KummerPoint {
    let product = p.iter().fold(Fp2::new(Fp::new([1, 0, 0, 0]), Fp::zero()), |acc, curr| acc.mul_mod(&curr, &PARAMS.p, PARAMS.mu));

    if !product.is_zero() {
        let pinv = p[3].inv_mod(&PARAMS.p, PARAMS.mu).unwrap();
        return [p[0].mul_mod(&pinv, &PARAMS.p, PARAMS.mu), p[1].mul_mod(&pinv, &PARAMS.p, PARAMS.mu), p[2].mul_mod(&pinv, &PARAMS.p, PARAMS.mu), Fp2::new(Fp::from_limb(1u64), Fp::zero())];
    } else {
        let pinv = p[0].inv_mod(&PARAMS.p, PARAMS.mu).unwrap();
        return [Fp2::new(Fp::from_limb(1u64), Fp::zero()), p[1].mul_mod(&pinv, &PARAMS.p, PARAMS.mu), p[2].mul_mod(&pinv, &PARAMS.p, PARAMS.mu), p[3].mul_mod(&pinv, &PARAMS.p, PARAMS.mu)];
    }
}

#[inline(always)]
pub fn pseudo_double(p: &KummerPoint, k: &[[Fp2; 4]; 4]) -> KummerPoint {
    four_way_mult(&hadamard(&four_way_mult(&squaring(&hadamard(&squaring(p))), &k[3])), &k[2])
}

#[inline(always)]
pub fn pseudo_add(p: &KummerPoint, q: &KummerPoint, r: &KummerPoint, k: &[[Fp2; 4]; 4]) -> KummerPoint {
    let hsp = hadamard(&squaring(p));
    let hsq = hadamard(&squaring(q));
    let s = four_way_mult(&hsp, &k[3]);
    four_way_mult(&hadamard(&four_way_mult(&hsq, &s)), &invert_4_constants(r))
}

#[inline(always)]
pub fn pseudo_triple(p: &KummerPoint, tc: &[[Fp2; 4]; 5]) -> KummerPoint {
    let r = hadamard(&squaring(p));
    let q = hadamard(&squaring(&four_way_mult(&hadamard(&four_way_mult(&squaring(&r), &tc[4])), &tc[3])));
    let s = four_way_mult(&r, &tc[4]);
    four_way_mult(&hadamard(&four_way_mult(&q, &s)), &invert_4_constants(p))
}

#[inline(always)]
pub const fn ind(i: &[[usize; 2]]) -> usize {
    // ind: 0b_0000_<I[1][1], I[1][0], I[0][1], I[0][0]>
    let l = (i[1][0] > i[0][0], i[1][1] > i[0][1]);
    match l {
        (false, false) => {
            0
        },
        (true, true) => {
            1
        },
        (true, false) => {
            2
        },
        (false, true) => {
            3
        }
    }
}

#[inline(always)]
pub fn double_thrice_add(r: &KummerPoint, s: &KummerPoint, del_rs: &KummerPoint, t: &KummerPoint, u: &KummerPoint, del_tu: &KummerPoint, v: &KummerPoint, w: &KummerPoint, del_vw: &KummerPoint, k: &[[Fp2; 4]; 4]) -> [KummerPoint; 4] {
    [pseudo_double(r, k), pseudo_add(r, s, del_rs, k), pseudo_add(t, u, del_tu, k), pseudo_add(v, w, del_vw, k)]
}

pub fn three_dac<const L: usize, Msg: UnitMessageEncodable>(p: &[[Fp2; 4]; 4], d: &[[Fp2; 4]; 4], dt: &[[Fp2; 4]; 4], m: &[Msg; L], n: &[Msg; L], k: &[[Fp2; 4]; 4]) -> [Fp2; 4] {
    let mut q = p.clone();

    let (first_add, [b0, b1, b2, b3]) = encode(m, n);

    let mut ii = [[1usize, 1usize], [2, 2], [1, 1]];

    if first_add {
        q[2] = pseudo_add(&q[2], &d[1], &d[0], k);
        ii[2][1] += 1;
    } else {
        q[2] = pseudo_add(&q[2], &d[0], &d[1], k);
        ii[2][0] += 1;
    }

    for i in (0..(Msg::SIZE * L - 1)).rev() {
        let ll = (b0[i], b1[i], b2[i], b3[i]);
        let ii0 = [ii[0][0] + ii[1][0], ii[0][1] + ii[1][1]];
        match ll {
            (false, false, false, false) => {
                ii = [
                    ii0,
                    [2 * ii[1][0], 2 * ii[1][1]],
                    [ii[1][0] + ii[2][0], ii[1][1] + ii[2][1]]
                ];
                
                [
                    q[1], q[0], q[2], q[3]
                ] = double_thrice_add(&q[1], &q[0], &d[2], &q[2], &q[1], &d[1], &q[1], &q[3], &dt[ind(&ii)], k);
            },
            (false, false, false, true) => {
                ii = [
                    ii0,
                    [2 * ii[1][0], 2 * ii[1][1]],
                    [ii[1][0] + ii[2][0], ii[1][1] + ii[2][1]]
                ];
                
                [
                    q[1], q[0], q[2], q[3]
                ] = double_thrice_add(&q[1], &q[0], &d[2], &q[2], &q[1], &d[0], &q[1], &q[3], &dt[ind(&ii)], k);
            },
            (false, false, true, false) => {
                ii = [
                    ii0,
                    [2 * ii[1][0], 2 * ii[1][1]],
                    [ii[1][0] + ii[2][0], ii[1][1] + ii[2][1]]
                ];
                
                [
                    q[1], q[0], q[2], q[3]
                ] = double_thrice_add(&q[1], &q[0], &d[3], &q[2], &q[1], &d[1], &q[1], &q[3], &dt[ind(&ii)], k);
            },
            (false, false, true, true) => {
                ii = [
                    ii0,
                    [2 * ii[1][0], 2 * ii[1][1]],
                    [ii[1][0] + ii[2][0], ii[1][1] + ii[2][1]]
                ];
                
                [
                    q[1], q[0], q[2], q[3]
                ] = double_thrice_add(&q[1], &q[0], &d[3], &q[2], &q[1], &d[0], &q[1], &q[3], &dt[ind(&ii)], k);
            },



            (false, true, false, false) => {
                ii = [
                    ii0,
                    [2 * ii[0][0], 2 * ii[0][1]],
                    [ii[0][0] + ii[2][0], ii[0][1] + ii[2][1]]
                ];
                
                [
                    q[1], q[0], q[2], q[3]
                ] = double_thrice_add(&q[0], &q[1], &d[2], &q[2], &q[0], &d[1], &q[0], &q[3], &dt[ind(&ii)], k);
            },
            (false, true, false, true) => {
                ii = [
                    ii0,
                    [2 * ii[0][0], 2 * ii[0][1]],
                    [ii[0][0] + ii[2][0], ii[0][1] + ii[2][1]]
                ];
                
                [
                    q[1], q[0], q[2], q[3]
                ] = double_thrice_add(&q[0], &q[1], &d[2], &q[2], &q[0], &d[0], &q[0], &q[3], &dt[ind(&ii)], k);
            },
            (false, true, true, false) => {
                ii = [
                    ii0,
                    [2 * ii[0][0], 2 * ii[0][1]],
                    [ii[0][0] + ii[2][0], ii[0][1] + ii[2][1]]
                ];
                
                [
                    q[1], q[0], q[2], q[3]
                ] = double_thrice_add(&q[0], &q[1], &d[3], &q[2], &q[0], &d[1], &q[0], &q[3], &dt[ind(&ii)], k);
            },
            (false, true, true, true) => {
                ii = [
                    ii0,
                    [2 * ii[0][0], 2 * ii[0][1]],
                    [ii[0][0] + ii[2][0], ii[0][1] + ii[2][1]]
                ];
                
                [
                    q[1], q[0], q[2], q[3]
                ] = double_thrice_add(&q[0], &q[1], &d[3], &q[2], &q[0], &d[0], &q[0], &q[3], &dt[ind(&ii)], k);
            },


            (true, false, false, false) => {
                ii = [
                    ii0,
                    [2 * ii[2][0], 2 * ii[2][1]],
                    [ii[1][0] + ii[2][0], ii[1][1] + ii[2][1]]
                ];
                [
                    q[1], q[2], q[0], q[3]
                ] = double_thrice_add(&q[2], &q[1], &d[1], &q[0], &q[1], &d[2], &q[2], &q[3], &dt[ind(&ii)], k);
            },
            (true, false, false, true) => {
                ii = [
                    ii0,
                    [2 * ii[2][0], 2 * ii[2][1]],
                    [ii[0][0] + ii[2][0], ii[0][1] + ii[2][1]]
                ];
                [
                    q[1], q[2], q[0], q[3]
                ] = double_thrice_add(&q[2], &q[0], &d[0], &q[0], &q[1], &d[2], &q[2], &q[3], &dt[ind(&ii)], k);
            },
            (true, false, true, false) => {
                ii = [
                    ii0,
                    [2 * ii[2][0], 2 * ii[2][1]],
                    [ii[1][0] + ii[2][0], ii[1][1] + ii[2][1]]
                ];
                [
                    q[1], q[2], q[0], q[3]
                ] = double_thrice_add(&q[2], &q[1], &d[1], &q[0], &q[1], &d[3], &q[2], &q[3], &dt[ind(&ii)], k);
            },
            (true, false, true, true) => {
                ii = [
                    ii0,
                    [2 * ii[2][0], 2 * ii[2][1]],
                    [ii[0][0] + ii[2][0], ii[0][1] + ii[2][1]]
                ];
                [
                    q[1], q[2], q[0], q[3]
                ] = double_thrice_add(&q[2], &q[0], &d[0], &q[0], &q[1], &d[3], &q[2], &q[3], &dt[ind(&ii)], k);
            },
            (true, true, false, false) => {
                ii = [
                    ii0,
                    [2 * ii[2][0], 2 * ii[2][1]],
                    [ii[0][0] + ii[2][0], ii[0][1] + ii[2][1]]
                ];
                [
                    q[1], q[2], q[0], q[3]
                ] = double_thrice_add(&q[2], &q[0], &d[1], &q[0], &q[1], &d[2], &q[2], &q[3], &dt[ind(&ii)], k);
            },
            (true, true, false, true) => {
                ii = [
                    ii0,
                    [2 * ii[2][0], 2 * ii[2][1]],
                    [ii[1][0] + ii[2][0], ii[1][1] + ii[2][1]]
                ];
                [
                    q[1], q[2], q[0], q[3]
                ] = double_thrice_add(&q[2], &q[1], &d[0], &q[0], &q[1], &d[2], &q[2], &q[3], &dt[ind(&ii)], k);
            },
            (true, true, true, false) => {
                ii = [
                    ii0,
                    [2 * ii[2][0], 2 * ii[2][1]],
                    [ii[0][0] + ii[2][0], ii[0][1] + ii[2][1]]
                ];
                [
                    q[1], q[2], q[0], q[3]
                ] = double_thrice_add(&q[2], &q[0], &d[1], &q[0], &q[1], &d[3], &q[2], &q[3], &dt[ind(&ii)], k);
            },
            (true, true, true, true) => {
                ii = [
                    ii0,
                    [2 * ii[2][0], 2 * ii[2][1]],
                    [ii[1][0] + ii[2][0], ii[1][1] + ii[2][1]]
                ];
                [
                    q[1], q[2], q[0], q[3]
                ] = double_thrice_add(&q[2], &q[1], &d[0], &q[0], &q[1], &d[3], &q[2], &q[3], &dt[ind(&ii)], k);
            }
        }
    }
    q[3]
}

