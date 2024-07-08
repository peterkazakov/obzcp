mod utils;

use std::time::Instant;

use rayon::prelude::*;
use std::cmp;
use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::io::{self, Write};
use std::sync::atomic::{AtomicUsize, Ordering};

// User must modify following two constants to get results for corresponding OBZCPs
const N: i16 = 29;
const MAX_ACC: i16 = 1;
// End of modification parameters

#[allow(dead_code)]
const P: i16 = 18;

const MASK: u64 = (1u64 << N) - 1;
const ZCZ: i16 = (N + 1) / 2;

const MIDDLE_ONE: u64 = 1u64 << (ZCZ - 1);
const MASK1: u64 = 0x2;
const MASK2: u64 = 1u64 << (N - 2);

const CHUNK_SIZE: i16 = (N - 3) / 2;
const NUMBER_OF_CHUNKS: u64 = 1u64 << CHUNK_SIZE;

#[allow(arithmetic_overflow)]
#[inline(always)]
fn conv_map(digit: u64) -> u64 {
    match digit {
        0x2 => 1u64 << (N - 2),
        0x4 => 1u64 << (N - 3),
        0x8 => 1u64 << (N - 4),
        0x10 => 1u64 << (N - 5),
        0x20 => 1u64 << (N - 6),
        0x40 => 1u64 << (N - 7),
        0x80 => 1u64 << (N - 8),
        0x100 => 1u64 << (N - 9),
        0x200 => 1u64 << (N - 10),
        0x400 => 1u64 << (N - 11),
        0x800 => 1u64 << (N - 12),
        0x1000 => 1u64 << (N - 13),
        0x2000 => 1u64 << (N - 14),
        0x4000 => 1u64 << (N - 15),
        0x8000 => 1u64 << (N - 16),
        0x10000 => 1u64 << (N - 17),
        0x20000 => 1u64 << (N - 18),
        0x40000 => 1u64 << (N - 19),
        0x80000 => 1u64 << (N - 20),
        0x100000 => 1u64 << (N - 21),
        0x200000 => 1u64 << (N - 22),
        0x400000 => 1u64 << (N - 23),
        0x800000 => 1u64 << (N - 24),
        0x1000000 => 1u64 << (N - 25),
        0x2000000 => 1u64 << (N - 26),
        _ => panic!("Unsupported value {:?}", digit),
    }
}

#[inline(always)]
fn get_weight(x: u64) -> i16 {
    x.count_ones() as i16
}

#[inline(always)]
fn hash(v: u64) -> u128 {
    let mut r = 0u128;
    let mut mask2: u64 = 0xffffffffffffffff;
    mask2 <<= 1;
    for t in 1..ZCZ {
        let c = get_weight(MASK & mask2 & (v ^ (v << t)));
        mask2 <<= 1;
        r = (r << 4) ^ (c as u128);
    }
    r
}

#[inline(always)]
fn reverse_hash(v: u64) -> u128 {
    let mut r = 0u128;
    let mut mask2: u64 = 0xffffffffffffffff;
    mask2 <<= 1;
    for t in 1..ZCZ {
        let c = get_weight(MASK & mask2 & (v ^ (v << t)));
        let ppb: i16 = N - t - c;
        mask2 <<= 1;
        r = (r << 4) ^ (ppb as u128);
    }
    r
}

#[inline(always)]
fn reverse(rr: u64, x1: u64) -> u64 {
    let mut res = 0u64;
    let mut m1 = 1u64;
    let mut m2 = MASK2;
    for _c in 0..((N - 2) / 2) {
        let x1 = (rr & m1) != (m1 & x1);
        if x1 {
            res |= m2;
        }
        m1 <<= 1;
        m2 >>= 1;
    }
    res
}

fn reverse_item(n: i16, x1: u64) -> u64 {
    let mut res = 0u64;
    let mut x = x1;
    for _c in 0..n {
        res <<= 1;
        res = res | (x & 1);
        x >>= 1;
    }
    res
}

#[allow(dead_code)]
fn inv(n: &i16, x1: &u64) -> u64 {
    x1 ^ ((1u64 << n) - 1)
}

#[allow(dead_code)]
fn rev(n: &i16, x1: &u64) -> u64 {
    reverse_item(*n, *x1)
}

#[allow(dead_code)]
fn is_solution(n: i16, max_acc: i16, a: u64, b: u64) -> bool {
    let mask: u64 = (1u64 << n) - 1;
    let zcz: i16 = (n + 1) / 2;
    let mut flag = true;

    let mut mask2: u64 = 0xffffffffffffffff;
    mask2 <<= 1;

    for t in 1..zcz {
        let ppb = get_weight(mask & mask2 & (a ^ (a << t)));
        let acc1: i16 = n - t - get_weight(mask & mask2 & (b ^ (b << t))) - ppb;
        if acc1 != 0 {
            flag = false;
            break;
        }
        mask2 <<= 1;
    }

    for t in zcz..n {
        let ppb = get_weight(mask & mask2 & (a ^ (a << t)));
        let acc1: i16 = n - t - get_weight(mask & mask2 & (b ^ (b << t))) - ppb;

        if acc1.abs() > max_acc {
            flag = false;
            break;
        }
        mask2 <<= 1;
    }
    flag
}

fn obzcp(mid_a: u64, mid_b: u64, odd: u64, pq_pairs: &HashMap<i16, Vec<i16>>) -> usize {
    let global_thread_count: AtomicUsize = AtomicUsize::new(0);
    let min_possible_weight: i16 = pq_pairs.keys().map(|v| *v).min().unwrap();
    let max_possible_weight: i16 = N - min_possible_weight;

    (0u64..NUMBER_OF_CHUNKS).into_par_iter().for_each(|x1| { //x1-->c
        let mut hash_codes: HashMap<u128, Vec<u64>> = HashMap::new();
        for c in 0..1u64 << CHUNK_SIZE {
            let a_sequence = (1u64 << (N - 1)) | mid_a | odd | (c << 1) | (reverse(c, x1));
            let p = get_weight(a_sequence);
            if (p >= min_possible_weight) && (p <= max_possible_weight) {
                let ro_a_seq = hash(a_sequence);
                match hash_codes.entry(ro_a_seq) {
                    Entry::Vacant(e) => { e.insert(vec![a_sequence]); },
                    Entry::Occupied(mut e) => { e.get_mut().push(a_sequence); },
                }
            }
        };

        let mut b_sequence: u64 = odd | mid_b;

        // Gray code counting
        let bit_to_change = 2u64; //shifted with one position
        b_sequence ^= bit_to_change;

        // Synchronize b's higher part with the lower's one
        let rr = (x1 << 1) ^ b_sequence;
        let mut m1 = MASK1;
        let mut m2 = MASK2;
        for _c in 1..=((N - 2) / 2) {
            let x1 = (rr & m1) != 0;
            let x2 = (rr & m2) != 0;
            if !(x1 ^ x2) {
                b_sequence |= m2;
            }
            m1 <<= 1;
            m2 >>= 1;
        }
        b_sequence = b_sequence | (1u64 << (N - 1));

        //iterate over all possible b sequenceswith Gray counting
        let number_iter:u64 = MIDDLE_ONE >> 1;
        for j in 2..number_iter
        {
            //Gray code
            let gray_idx: i64 = j as i64;
            let original_bit_to_change = (gray_idx & (-gray_idx)) as u64;
            let bit_to_change = original_bit_to_change << 1; //shifted with one position

            //update only 1 bit
            let conv_bit_to_change = conv_map(bit_to_change);
            b_sequence ^= conv_bit_to_change ^ bit_to_change;

            let a_values = hash_codes.get(&reverse_hash(b_sequence));
            if a_values != None {
                for a_sequence in a_values.unwrap().into_iter()
                {
                    if *a_sequence > b_sequence { //some ordering to skip duplication
                        let mut max = 2;
                        let mut mask2: u64 = 0xffffffffffffffff;
                        mask2 <<= ZCZ;
                        let mut max_acc= 2;
                        for t in ZCZ..N
                        {
                            let ppb = get_weight(MASK & mask2 & (b_sequence ^ (b_sequence << t)));
                            let acc1: i16 = N - t - get_weight(MASK & mask2 & (*a_sequence ^ (*a_sequence << t))) - ppb;
                            mask2 <<= 1;

                            if acc1.abs() > MAX_ACC {
                                max = 30;
                                break;
                            }
                            else {
                                max_acc = cmp::max(max_acc, 2*acc1.abs());
                            }
                        }

                        if max < 30 { //check the optimality
                            println!("['{:X}','{:X}', '{:#b}', '{:#b}']", a_sequence, b_sequence, a_sequence, b_sequence);
                            // do verification of the upper part
                            let mut mask2:u64 = 0xffffffffffffffff;
                            mask2 <<= 1;
                            // Be paranoic and recheck because of hashes (ro's in the publication)
                            for t in 1 .. ZCZ
                            {
                                let ppb = get_weight(MASK & mask2 & (b_sequence ^ (b_sequence << t)));
                                let acc1: i16 = N - t - get_weight(MASK & mask2 & (*a_sequence ^ (*a_sequence << t))) - ppb;
                                if acc1 != 0 {
                                    // Informative, indicates quality of the hash code construction
                                    println!("Oooooops, hash breaks on t={:}, acc={:}, hash={:#x}, reserve={:#x}!", t, acc1, hash(*a_sequence), reverse_hash(b_sequence));
                                }
                                mask2 <<= 1;
                            }

                            io::stdout().flush().unwrap();
                            global_thread_count.fetch_add(1, Ordering::SeqCst);
                        }
                    }
                }
            }
        }//a_values
    }); // x1
    global_thread_count.load(Ordering::SeqCst)
}

#[allow(dead_code)]
fn get_hexagonal_index(x: u64) -> u64 {
    let mut i = 1u64;
    let mut r: u64 = 1;
    while r < x {
        i += 1;
        r = 2 * i * (2 * i - 1) / 2;
    }
    if r == x {
        i
    } else {
        0
    }
}

#[allow(dead_code)]
fn is_equivalent(n: i16, a1: u64, b1: u64, a2: u64, b2: u64) -> bool {
    let x1 = match a1 {
        a if a == b2 => true,
        a if a == a2 => true,
        a if a == inv(&n, &b2) => true,
        a if a == rev(&n, &inv(&n, &b2)) => true,
        a if a == rev(&n, &b2) => true,
        a if a == inv(&n, &a2) => true,
        a if a == rev(&n, &inv(&n, &a2)) => true,
        a if a == rev(&n, &a2) => true,
        _ => false,
    };
    let x2 = match b1 {
        b if b == b2 => true,
        b if b == a2 => true,
        b if b == inv(&n, &b2) => true,
        b if b == rev(&n, &inv(&n, &b2)) => true,
        b if b == rev(&n, &b2) => true,
        b if b == inv(&n, &a2) => true,
        b if b == rev(&n, &inv(&n, &a2)) => true,
        b if b == rev(&n, &a2) => true,
        _ => false,
    };
    x1 && x2
}

#[allow(dead_code)]
fn is_h_regular(n: i16, a1: u64, b1: u64) -> bool {
    let mut res = true;
    let mut mask2: u64 = 1u64 << (n - 1);
    let mut mask1 = 1u64;
    let size = (n as i16 - 3) / 2;

    let c = a1 ^ b1;
    for _ in 1..=size {
        mask1 <<= 1;
        mask2 >>= 1;
        if (c & mask1 == 0) ^ (c & mask2 == 0) == true {
            res = false;
            break;
        }
    }
    res
}

#[allow(dead_code)]
fn is_h_reversed_regular(n: i16, a1: u64, b1: u64) -> bool {
    let mut res = true;
    let size = (n as i16 - 3) / 2;
    let mut mask2: u64 = 1u64 << (size + 1);
    let mut mask1 = 1u64;

    let c = a1 ^ b1;
    for _ in 1..=size {
        mask1 <<= 1;
        mask2 <<= 1;
        if (c & mask1 == 0) ^ (c & mask2 == 0) == false {
            res = false;
            break;
        }
    }
    res
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sol_25() {
        assert_eq!(
            is_solution(25, 1, reverse_item(25, 0x159fe24) ^ 0x1FFFFFF, 0x148e0ca),
            true
        );
    }

    #[test]
    fn test_sol_35() {
        assert_eq!(is_solution(35, 3, 0x7905A9444, 0x710C1A3B2), true);
    }

    #[test]
    fn test_sol_49() {
        assert_eq!(
            is_solution(
                49,
                1,
                0b1010110111100110100011111101011110011001101000000,
                0b1000001011001100111101010000011101001100001001010
            ),
            true
        );
    }

    #[test]
    fn test_h_reserved_49() {
        assert_eq!(
            is_h_reversed_regular(49, 0x156CC7FF8E494, 0x1524E3E03992A),
            true
        );
    }
}

fn main() {
    println!(
        "N={:}, max_acc={:}, chunks={:}",
        N,
        MAX_ACC * 2,
        NUMBER_OF_CHUNKS
    );

    rayon::ThreadPoolBuilder::new()
        .num_threads(16)
        .build_global()
        .unwrap();
    let pq_pairs = utils::calc_pq(N as i16);
    let start = Instant::now();

    let x = obzcp(MIDDLE_ONE, MIDDLE_ONE, 0, &pq_pairs);
    let y = obzcp(MIDDLE_ONE, MIDDLE_ONE, 1, &pq_pairs);
    let x2 = obzcp(0, MIDDLE_ONE, 0, &pq_pairs);
    let y2 = obzcp(0, MIDDLE_ONE, 1, &pq_pairs);
    let x3 = obzcp(MIDDLE_ONE, 0, 0, &pq_pairs);
    let y3 = obzcp(MIDDLE_ONE, 0, 1, &pq_pairs);
    let x4 = obzcp(0, 0, 0, &pq_pairs);
    let y4 = obzcp(0, 0, 1, &pq_pairs);

    let duration = start.elapsed();
    println!(
        "Found {:?} solutions for {:?}.",
        x + y + x2 + y2 + x3 + y3 + x4 + y4,
        duration
    );
}
