use crate::field_types::*;
use crate::isogeny::get_isogeny_33_image_and_points;
use crate::kummer::get_triple_constants_from_fund_thetas;
use crate::kummer::pseudo_triple;
use crate::kummer::three_dac;
use crate::params::PARAMS;
use crate::kummer;

pub struct KummerHashSecp256;

impl KummerHashSecp256 {
    pub fn new() -> Self { Self }

    // hash_optimal 포팅
    pub fn hash_optimal(&self, msg: &[Fp; 3]) -> [Fp2; 4] {
        let mut tc = kummer::get_triple_constants_from_fund_thetas(&PARAMS.k[1]);

        let qq1 = three_dac(&PARAMS.gens[0], &PARAMS.gens[1], &PARAMS.gens[2], &msg[0].get_limbs(), &msg[1].get_limbs(), &PARAMS.k);
        let qq2 = three_dac(&PARAMS.gens[3], &PARAMS.gens[4], &PARAMS.gens[5], &msg[1].get_limbs(), &msg[2].get_limbs(), &PARAMS.k);

        let mut pts = Vec::<(KummerPoint, KummerPoint)>::new();
        let mut ind = 0usize;
        let mut inds = Vec::<usize>::new();
        let mut image_thetas = [Fp2::zero(); 4];

        let (mut rr, mut ss) = (qq1, qq2);

        for row in 1..PARAMS.e {
            while ind < (PARAMS.e - row) as usize {
                pts.push((rr, ss));
                inds.push(ind);

                let m = PARAMS.strategy[PARAMS.e as usize - ind - row as usize];

                for _ in 0..m {
                    (rr, ss) = (pseudo_triple(&rr, &tc), pseudo_triple(&ss, &tc));
                }
                ind += m;
            }

            (pts, image_thetas) = get_isogeny_33_image_and_points(&rr, &ss, &tc, &pts);
            tc = get_triple_constants_from_fund_thetas(&image_thetas);

            if pts.len() > 0 {
                assert_eq!(inds.len(), pts.len());
                (rr, ss) = pts.pop().unwrap();
                ind = inds.pop().unwrap();
            }
        }

        image_thetas
    }
}