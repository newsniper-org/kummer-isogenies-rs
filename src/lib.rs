#![feature(generic_const_exprs)]
#![feature(const_trait_impl)]
#![feature(slice_as_array)]

use num_bigint::BigUint;
use num_traits::FromBytes;

use crate::field_types::Fp;

pub fn biguint_to_fp(val: BigUint) -> Fp {
    let digits = val.to_u64_digits();
    let mut result = [0u64; 4];
    let minlen = digits.len().min(4usize);
    for i in 0..minlen {
        result[i] = digits[i];
    }
    Fp::new(result)    
}

pub fn fp_to_biguint(val: Fp) -> BigUint {
    let mut le_bytes = [0u8; 32];
    val.to_le_bytes(&mut le_bytes);
    BigUint::from_le_bytes(&le_bytes)
}

pub mod field_types;
pub mod params;
pub mod kummer;

pub mod coding;
pub mod isogeny;
pub mod kuhash;