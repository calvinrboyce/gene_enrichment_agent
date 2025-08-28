import numpy as np
import torch
from tqdm.auto import tqdm
from transformers import AutoTokenizer, AutoModel

import warnings
import os
warnings.filterwarnings("ignore", category=FutureWarning, module="torch.nn.modules.module")
os.environ["TOKENIZERS_PARALLELISM"] = "false"

# get device
if torch.backends.mps.is_available():
    device = "mps"
elif torch.cuda.is_available():
    device = "cuda"
else:
    device = "cpu"

# load tokenizer and model
tokenizer = AutoTokenizer.from_pretrained("cambridgeltl/SapBERT-from-PubMedBERT-fulltext")  
model = AutoModel.from_pretrained("cambridgeltl/SapBERT-from-PubMedBERT-fulltext").to(device)

def embed_terms(terms):
    """
    Embed terms using SapBERT.
    """

    # embed terms
    bs = 32 # batch size during inference
    all_embs = []
    for i in np.arange(0, len(terms), bs):
        toks = tokenizer.batch_encode_plus(terms[i:i+bs], 
                                        padding="max_length", 
                                        max_length=25, 
                                        truncation=True,
                                        return_tensors="pt")
        toks_cuda = {}
        for k,v in toks.items():
            toks_cuda[k] = v.to(device)
        cls_rep = model(**toks_cuda)[0][:,0,:] # use CLS representation as the embedding
        all_embs.append(cls_rep.cpu().detach().numpy())

    all_embs = np.concatenate(all_embs, axis=0)

    return all_embs
