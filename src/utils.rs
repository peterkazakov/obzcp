
use std::cmp;
use std::collections::{HashMap, HashSet};

pub fn calc_pq(n: i16) -> HashMap<i16, Vec<i16>> {
    let mut d = 1i16;
    let mut f = 0;

    let n2 = ((2 * n - 1) as f64).sqrt() as i16;
    let mut pq = vec![vec![0i16; 250]; 250];
    pq[0][0] = 0;
    pq[0][1] = 0;
    for i in 1..=n {
        for j in 1..=n {
            let a = (i - j).abs();
            let b = (n - i - j).abs();
            if (a <= n2) && (b <= n2) && ((a > 0) || (b > 0)) {
                pq[d as usize][d as usize] = i;
                pq[d as usize][(d + 1) as usize] = j;
                d = d + 1;
                f = d;
            }
        }
    }

    let mut pq2 = HashSet::new();
    for i in 0..f {
        let a2 = cmp::min(pq[i as usize][i as usize], n - pq[i as usize][i as usize]);
        let b2 = cmp::min(
            pq[i as usize][(i + 1) as usize],
            n - pq[i as usize][(i + 1) as usize],
        );

        if a2 <= b2 {
            pq2.insert((a2, b2));
        }
    }

    let mut pqfix = HashMap::new();
    for (k, v) in pq2 {
        pqfix.entry(k).or_insert_with(Vec::new).push(v)
    }
    println!(" {:?}", pqfix);
    pqfix
}
