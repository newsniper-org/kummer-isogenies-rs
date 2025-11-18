pub trait MessageEncodable {
    fn encode(&self) -> Vec<bool>;
}

pub trait UnitMessageEncodable : Sized + MessageEncodable {
    const SIZE: usize;
}

impl MessageEncodable for u8 {
    fn encode(&self) -> Vec<bool> {
        [
            0b00000001u8, 0b00000010, 0b00000100, 0b00001000,0b00010000, 0b00100000, 0b01000000, 0b10000000
        ].map(|flag|(self & flag) != 0u8).to_vec()
    }
}
impl UnitMessageEncodable for u8 {
    const SIZE: usize = 8;
}

impl<T: MessageEncodable> MessageEncodable for [T] {
    fn encode(&self) -> Vec<bool> {
        self.iter().map(|unit|unit.encode()).collect::<Vec<Vec<bool>>>().concat()
    }
}


pub fn encode<const L: usize, Msg: UnitMessageEncodable>(m: &[Msg; L], n: &[Msg; L]) -> (bool, [Vec<bool>; 4]) {
    let m_encoded = m.encode();
    let n_encoded = n.encode();

    let (
        mut b0,
        mut b1,
        mut b2,
        mut b3
    ) = (
        Vec::<bool>::new(),
        Vec::<bool>::new(),
        Vec::<bool>::new(),
        Vec::<bool>::new()
    );

    let mut first_add = m_encoded[0];

    for i in 1..(Msg::SIZE*L) {
        let mm = m_encoded[i-1] ^ m_encoded[i];
        let nn = n_encoded[i-1] ^ n_encoded[i];
        let mn = mm ^ nn;
        b0.push(mn);
        b1.push(mm);
        b2.push(m_encoded[i] ^ n_encoded[i]);
        b3.push(first_add);
        first_add = mm ^ (!mn && first_add);
    }
    (first_add, [b0, b1, b2, b3])
}


impl MessageEncodable for u64 {
    fn encode(&self) -> Vec<bool> {
        let bytes: [u8; 8] = self.to_le_bytes();
        bytes.encode()
    }
}

impl UnitMessageEncodable for u64 {
    const SIZE: usize = 64;
}