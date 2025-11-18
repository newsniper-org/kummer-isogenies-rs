pub const N_LIMBS: usize = 4; // 4 * 64 = 256 bits

pub const EXP: usize = 154;

/// F_p 원소 (252비트 소수)
pub type Fp = cubemoma::BigField<N_LIMBS>;

/// F_p^2 원소 (i^2 = -1 기반)
/// cubemoma::FpComplex가 이 구조와 일치
pub type Fp2 = cubemoma::FpComplex<N_LIMBS>;

/// 쿠머 곡면 위의 점 (프로젝티브 좌표)
pub type KummerPoint = [Fp2; 4];