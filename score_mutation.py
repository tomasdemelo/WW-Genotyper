from transformers import AutoTokenizer, AutoModelForMaskedLM
import torch

_tokenizer = None
_model = None

def get_esm_model():
    global _tokenizer, _model
    if _tokenizer is None or _model is None:
        model_name = "facebook/esm2_t6_8M_UR50D"
        try:
            _tokenizer = AutoTokenizer.from_pretrained(model_name)
            _model = AutoModelForMaskedLM.from_pretrained(model_name)
        except Exception as e:
            return None, None
    return _tokenizer, _model

def calculate_mutation_risk(wt_sequence, mutated_position_1_based, mutant_aa):
    tokenizer, model = get_esm_model()
    if tokenizer is None or model is None:
        return 0.0 # Fail safe if offline or huggingface is down

    # Clean sequence just in case
    wt_sequence = wt_sequence.replace("*", "")
    
    if mutated_position_1_based > len(wt_sequence) or mutated_position_1_based < 1:
        return 0.0 # Invalid position

    seq_list = list(wt_sequence)
    wt_aa = seq_list[mutated_position_1_based - 1]
    
    # Mask the target position
    seq_list[mutated_position_1_based - 1] = tokenizer.mask_token
    masked_seq = "".join(seq_list)
    
    inputs = tokenizer(masked_seq, return_tensors="pt")
    
    with torch.no_grad():
        logits = model(**inputs).logits
        
    mask_token_index = (inputs.input_ids == tokenizer.mask_token_id)[0].nonzero(as_tuple=True)[0]
    
    if len(mask_token_index) == 0:
        return 0.0

    wt_token_id = tokenizer.convert_tokens_to_ids(wt_aa)
    mut_token_id = tokenizer.convert_tokens_to_ids(mutant_aa)
    
    wt_logit = logits[0, mask_token_index, wt_token_id].item()
    mut_logit = logits[0, mask_token_index, mut_token_id].item()
    
    # Negative log-odds means the mutation breaks evolutionary rules (high risk)
    log_odds = mut_logit - wt_logit
    
    return float(log_odds)
