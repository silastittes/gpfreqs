
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_p_vec() {
        let zero_probs = "0,0,0".to_string();
        p_vec(&zero_probs);
    }

    #[test]
    fn test_add() {
        assert_eq!(add(1, 2), 3);
    }

    #[test]
    fn test_bad_add() {
        // This assert would fire and test will fail.
        // Please note, that private functions can be tested too!
        assert_eq!(bad_add(1, 2), 3);
    }
}
