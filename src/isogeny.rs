use crate::{field_types::{Fp, Fp2, KummerPoint}, kummer::{four_way_mult, hadamard, squaring}, params::PARAMS};

pub fn isogeny_33_eval(p: &KummerPoint, coeffs: &[Fp2; 5]) -> [Fp2; 4] {
    let p_sq = squaring(p);

    let u = coeffs;
    let [x, y, z, t] = *p;

    let (u4xy, u4zt) = (
        u[4].mul_mod(&x, &PARAMS.p, PARAMS.mu).mul_mod(&y, &PARAMS.p, PARAMS.mu),
        u[4].mul_mod(&z, &PARAMS.p, PARAMS.mu).mul_mod(&t, &PARAMS.p, PARAMS.mu)
    );

    let p1 = four_way_mult(&[u[0], u[1], u[2], u[3]], &p_sq)
        .iter().fold(Fp2::zero(), |acc, curr| acc.add_mod(&curr, &PARAMS.p))
        .mul_mod(&p[0], &PARAMS.p, PARAMS.mu)
        .add_mod(&u4zt.mul_mod(&p[1], &PARAMS.p, PARAMS.mu), &PARAMS.p);

    let p2 = four_way_mult(&[u[0], u[1], u[2], u[3]], &[p_sq[1], p_sq[0], p_sq[3], p_sq[2]])
        .iter().fold(Fp2::zero(), |acc, curr| acc.add_mod(&curr, &PARAMS.p))
        .mul_mod(&p[1], &PARAMS.p, PARAMS.mu)
        .add_mod(&u4zt.mul_mod(&p[0], &PARAMS.p, PARAMS.mu), &PARAMS.p);

    let p3 = four_way_mult(&[u[0], u[1], u[2], u[3]], &[p_sq[2], p_sq[3], p_sq[0], p_sq[1]])
        .iter().fold(Fp2::zero(), |acc, curr| acc.add_mod(&curr, &PARAMS.p))
        .mul_mod(&p[2], &PARAMS.p, PARAMS.mu)
        .add_mod(&u4xy.mul_mod(&p[3], &PARAMS.p, PARAMS.mu), &PARAMS.p);

    let p4 = four_way_mult(&[u[0], u[1], u[2], u[3]], &[p_sq[3], p_sq[2], p_sq[1], p_sq[0]])
        .iter().fold(Fp2::zero(), |acc, curr| acc.add_mod(&curr, &PARAMS.p))
        .mul_mod(&p[3], &PARAMS.p, PARAMS.mu)
        .add_mod(&u4xy.mul_mod(&p[2], &PARAMS.p, PARAMS.mu), &PARAMS.p);

    [p1, p2, p3, p4]
}


pub fn get_isogeny_33_image_and_points(r: &KummerPoint, s: &KummerPoint, tc: &[[Fp2; 4]; 5], pts: &[(KummerPoint, KummerPoint)]) -> (Vec<([Fp2; 4], [Fp2; 4])>, [Fp2; 4]) {
    let [a, b, c, d] = tc[0];
    let [xr, yr, zr, tr] = *r;
    let [xs, ys, zs, ts] = *s;
    let s0 = tc[1];
    let h0_inv = tc[4];

    let (ab, ac, ad, bc, bd, cd) = (
        a.mul_mod(&b, &PARAMS.p, PARAMS.mu),
        a.mul_mod(&c, &PARAMS.p, PARAMS.mu),
        a.mul_mod(&d, &PARAMS.p, PARAMS.mu),
        b.mul_mod(&c, &PARAMS.p, PARAMS.mu),
        b.mul_mod(&d, &PARAMS.p, PARAMS.mu),
        c.mul_mod(&d, &PARAMS.p, PARAMS.mu)
    );

    let (
        dd1p, dd1m,
        dd2p, dd2m,
        dd3p, dd3m
    ) = (
        ab.add_mod(&cd, &PARAMS.p),
        ab.sub_mod(&cd, &PARAMS.p),
        ac.add_mod(&bd, &PARAMS.p),
        ac.sub_mod(&bd, &PARAMS.p),
        ad.add_mod(&bc, &PARAMS.p),
        ad.sub_mod(&bc, &PARAMS.p)
    );

    let (dd1, dd2, dd3) = (
        dd1p.mul_mod(&dd1m, &PARAMS.p, PARAMS.mu),
        dd2p.mul_mod(&dd2m, &PARAMS.p, PARAMS.mu),
        dd3p.mul_mod(&dd3m, &PARAMS.p, PARAMS.mu)
    );

    let (dd12, dd13, dd23) =(
        
        dd1.mul_mod(&dd2, &PARAMS.p, PARAMS.mu),
        dd1.mul_mod(&dd3, &PARAMS.p, PARAMS.mu),
        dd2.mul_mod(&dd3, &PARAMS.p, PARAMS.mu)
    );

    let (hsr, hss) = (
        hadamard(&squaring(r)),
        hadamard(&squaring(s))
    );

    let (hr0, hs0) = (
        four_way_mult(&hsr, &h0_inv),
        four_way_mult(&hss, &h0_inv)
    );

    let (hr, hs) = (
        hadamard(&hr0), hadamard(&hs0)
    );

    let [rxx, ryy, rzz, rtt] = hr;
    let [sxx, syy, szz, stt] = hs;

    let (xyr, xzr, xtr, yzr, ytr, ztr) = (
        xr.mul_mod(&yr, &PARAMS.p, PARAMS.mu),
        xr.mul_mod(&zr, &PARAMS.p, PARAMS.mu),
        xr.mul_mod(&tr, &PARAMS.p, PARAMS.mu),
        yr.mul_mod(&zr, &PARAMS.p, PARAMS.mu),
        yr.mul_mod(&tr, &PARAMS.p, PARAMS.mu),
        zr.mul_mod(&tr, &PARAMS.p, PARAMS.mu)
    );

    let (xys, xzs, xts, yzs, yts, zts) = (
        xs.mul_mod(&ys, &PARAMS.p, PARAMS.mu),
        xs.mul_mod(&zs, &PARAMS.p, PARAMS.mu),
        xs.mul_mod(&ts, &PARAMS.p, PARAMS.mu),
        ys.mul_mod(&zs, &PARAMS.p, PARAMS.mu),
        ys.mul_mod(&ts, &PARAMS.p, PARAMS.mu),
        zs.mul_mod(&ts, &PARAMS.p, PARAMS.mu)
    );

    
    let (mut rsig, mut rdel, mut ssig, mut sdel) = (
        xyr.add_mod(&ztr, &PARAMS.p).mul_mod(&dd1m, &PARAMS.p, PARAMS.mu),
        xyr.sub_mod(&ztr, &PARAMS.p).mul_mod(&dd1p, &PARAMS.p, PARAMS.mu),
        xys.add_mod(&zts, &PARAMS.p).mul_mod(&dd1m, &PARAMS.p, PARAMS.mu),
        xys.sub_mod(&zts, &PARAMS.p).mul_mod(&dd1p, &PARAMS.p, PARAMS.mu),
    );
    let (rxy, rzt, sxy, szt) = (
        rsig.add_mod(&rdel, &PARAMS.p),
        rsig.sub_mod(&rdel, &PARAMS.p),
        ssig.add_mod(&sdel, &PARAMS.p),
        ssig.sub_mod(&sdel, &PARAMS.p)
    );
    (rsig, rdel, ssig, sdel) = (
        xzr.add_mod(&ytr, &PARAMS.p).mul_mod(&dd2m, &PARAMS.p, PARAMS.mu),
        xzr.sub_mod(&ytr, &PARAMS.p).mul_mod(&dd2p, &PARAMS.p, PARAMS.mu),
        xzs.add_mod(&yts, &PARAMS.p).mul_mod(&dd2m, &PARAMS.p, PARAMS.mu),
        xzs.sub_mod(&yts, &PARAMS.p).mul_mod(&dd2p, &PARAMS.p, PARAMS.mu),
    );
    let (rxz, ryt, sxz, syt) = (
        rsig.add_mod(&rdel, &PARAMS.p),
        rsig.sub_mod(&rdel, &PARAMS.p),
        ssig.add_mod(&sdel, &PARAMS.p),
        ssig.sub_mod(&sdel, &PARAMS.p)
    );
    (rsig, rdel, ssig, sdel) = (
        xtr.add_mod(&yzr, &PARAMS.p).mul_mod(&dd3m, &PARAMS.p, PARAMS.mu),
        xtr.sub_mod(&yzr, &PARAMS.p).mul_mod(&dd3p, &PARAMS.p, PARAMS.mu),
        xts.add_mod(&yzs, &PARAMS.p).mul_mod(&dd3m, &PARAMS.p, PARAMS.mu),
        xts.sub_mod(&yzs, &PARAMS.p).mul_mod(&dd3p, &PARAMS.p, PARAMS.mu),
    );
    let (rxt, ryz, sxt, syz) = (
        rsig.add_mod(&rdel, &PARAMS.p),
        rsig.sub_mod(&rdel, &PARAMS.p),
        ssig.add_mod(&sdel, &PARAMS.p),
        ssig.sub_mod(&sdel, &PARAMS.p)
    );

    let c2 = {
        let tmp1 = [&dd23, &dd13, &dd12];
        let tmp2 = [&sxy, &sxz, &sxt];
        let result = [0usize, 1, 2].map(|idx| tmp1[idx].mul_mod(tmp2[idx], &PARAMS.p, PARAMS.mu)).iter().fold(Fp2::zero(), |acc, curr| acc.add_mod(&curr, &PARAMS.p));
        result
    };
    let d2 = {
        let tmp1 = [&dd23, &dd13, &dd12];
        let tmp2 = [&rxy, &rxz, &rxt];
        let result = [0usize, 1, 2].map(|idx| tmp1[idx].mul_mod(tmp2[idx], &PARAMS.p, PARAMS.mu)).iter().fold(Fp2::zero(), |acc, curr| acc.add_mod(&curr, &PARAMS.p));
        result
    };

    let ee1 = {
        let left = d2.mul_mod(&szt,&PARAMS.p, PARAMS.mu);
        let right = c2.mul_mod(&rzt,&PARAMS.p, PARAMS.mu);
        left.sub_mod(&right, &PARAMS.p).mul_mod(&dd23, &PARAMS.p, PARAMS.mu)
    };

    let ee2 = {
        let left = sxx.mul_mod(&ryy,&PARAMS.p, PARAMS.mu);
        let right = rxx.mul_mod(&syy,&PARAMS.p, PARAMS.mu);
        left.sub_mod(&right, &PARAMS.p)
    };

    let (ee1r, ee1s, ee2c, ee2d) = (
        ee1.mul_mod(&rxx, &PARAMS.p, PARAMS.mu),
        ee1.mul_mod(&sxx, &PARAMS.p, PARAMS.mu),
        ee2.mul_mod(&c2, &PARAMS.p, PARAMS.mu),
        ee2.mul_mod(&d2, &PARAMS.p, PARAMS.mu)
    );

    let two = Fp2::new(Fp::from_limb(2u64), Fp::zero());

    let u1 = ee1r.mul_mod(&sxx, &PARAMS.p, PARAMS.mu).mul_mod(&two,  &PARAMS.p, PARAMS.mu);
    let u2 = {
        let interm = ee1s.mul_mod(&ryy, &PARAMS.p, PARAMS.mu).add_mod(&ee1r.mul_mod(&syy, &PARAMS.p, PARAMS.mu), &PARAMS.p);
        let tmp = ee2d.mul_mod(&szt, &PARAMS.p, PARAMS.mu).add_mod(&ee2c.mul_mod(&rzt, &PARAMS.p, PARAMS.mu), &PARAMS.p).mul_mod(&dd23, &PARAMS.p, PARAMS.mu);
        interm.add_mod(&tmp, &PARAMS.p)
    };
    let u3 = {
        let interm = ee1s.mul_mod(&rzz, &PARAMS.p, PARAMS.mu).add_mod(&ee1r.mul_mod(&szz, &PARAMS.p, PARAMS.mu), &PARAMS.p);
        let tmp = ee2d.mul_mod(&syt, &PARAMS.p, PARAMS.mu).add_mod(&ee2c.mul_mod(&ryt, &PARAMS.p, PARAMS.mu), &PARAMS.p).mul_mod(&dd13, &PARAMS.p, PARAMS.mu);
        interm.add_mod(&tmp, &PARAMS.p)
    };
    let u4 = {
        let interm = ee1s.mul_mod(&rtt, &PARAMS.p, PARAMS.mu).add_mod(&ee1r.mul_mod(&stt, &PARAMS.p, PARAMS.mu), &PARAMS.p);
        let tmp = ee2d.mul_mod(&syz, &PARAMS.p, PARAMS.mu).add_mod(&ee2c.mul_mod(&ryz, &PARAMS.p, PARAMS.mu), &PARAMS.p).mul_mod(&dd12, &PARAMS.p, PARAMS.mu);
        interm.add_mod(&tmp, &PARAMS.p)
    };
    let u5 = ee2c.mul_mod(&d2, &PARAMS.p, PARAMS.mu).mul_mod(&two, &PARAMS.p, PARAMS.mu);

    let coeffs = [
        u1, u2, u3, u4, u5
    ];

    let (u5ab, u5cd) = (
        u5.mul_mod(&ab, &PARAMS.p, PARAMS.mu),
        u5.mul_mod(&cd, &PARAMS.p, PARAMS.mu)
    );

    let [s0a2, s0b2, s0c2, s0d2] = s0;

    let aa = {
        let tmp = four_way_mult(&[u1, u2, u3, u4], &[s0a2, s0b2, s0c2, s0d2]).iter().fold(Fp2::zero(), |acc, curr| acc.sub_mod(&curr, &PARAMS.p)).mul_mod(&a, &PARAMS.p, PARAMS.mu);
        u5cd.mul_mod(&b, &PARAMS.p, PARAMS.mu).add_mod(&tmp, &PARAMS.p)
    };
    let bb = {
        let tmp = four_way_mult(&[u1, u2, u3, u4], &[s0b2, s0a2, s0d2, s0c2]).iter().fold(Fp2::zero(), |acc, curr| acc.sub_mod(&curr, &PARAMS.p)).mul_mod(&a, &PARAMS.p, PARAMS.mu);
        u5cd.mul_mod(&a, &PARAMS.p, PARAMS.mu).add_mod(&tmp, &PARAMS.p)
    };
    let cc = {
        let tmp = four_way_mult(&[u1, u2, u3, u4], &[s0c2, s0d2, s0a2, s0b2]).iter().fold(Fp2::zero(), |acc, curr| acc.sub_mod(&curr, &PARAMS.p)).mul_mod(&a, &PARAMS.p, PARAMS.mu);
        u5ab.mul_mod(&d, &PARAMS.p, PARAMS.mu).add_mod(&tmp, &PARAMS.p)
    };
    let dd = {
        let tmp = four_way_mult(&[u1, u2, u3, u4], &[s0d2, s0c2, s0b2, s0a2]).iter().fold(Fp2::zero(), |acc, curr| acc.sub_mod(&curr, &PARAMS.p)).mul_mod(&a, &PARAMS.p, PARAMS.mu);
        u5ab.mul_mod(&c, &PARAMS.p, PARAMS.mu).add_mod(&tmp, &PARAMS.p)
    };

    let image_thetas = [aa, bb, cc, dd];

    let pushed_pts = pts.iter().map(|(pt0, pt1)| {
        (isogeny_33_eval(&pt0, &coeffs), isogeny_33_eval(&pt1, &coeffs))
    }).collect::<Vec<(KummerPoint, KummerPoint)>>();

    (pushed_pts, image_thetas)
}