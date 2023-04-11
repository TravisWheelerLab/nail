use nale::align::needleman_wunsch::{needleman_wunsch, print_alignment, print_trace};
use nale::structs::Sequence;
use std::str::FromStr;

pub fn main() {
    let s1 = Sequence::from_utf8(
        "SGGMDIEGAFAKVDRSVSTGSVPPPAPEDTPTLSDQAAVRNAVSSADTKQLKDSPVPWANPTTGSAGVITPIAEAYANGRTCRQFLTSRHAFDGIAWFQGRTCRLGQGQWQLTSFRPW"
            .as_bytes(),
    )
    .unwrap();

    let s2 = Sequence::from_utf8(
        "gldiegavkvdrsvvtGsvaesaeededkesDqaavrnavssvdlkklkdssvPWenattGsaGvitsiaedsangklCRqFltsreafegiakfqGetClleqGeWeltslkpa"
            .as_bytes(),
    )
    .unwrap();

    let trace = needleman_wunsch(&s1, &s2);

    print_trace(&trace);
    print_alignment(&trace, &s1, &s2);
}
