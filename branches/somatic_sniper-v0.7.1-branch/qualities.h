#ifndef QUALITIES

//this is for creating a gigantic lookup table of mapping qualities and base qualities given a base

//according to the samtools spec base qualities range from 0 to 93 and mapping quality ranges from 0 to 254 with 255 meaning that mapping quality is not available

//for each base we have the probability P(o|b) 

//will be done in decimal to start

double probability_call_correct[4][4][94][256];

//TODO move to log space
////TODO split into properly compartmentalized files or even a hardcoded table
void initialize_call_probabilities() {
    int base = called_base = base_quality = map_quality = 0;
    for(base = 0; base < 4; base++) {
        for(called_base = 0; base < 4; base++) {
            for(base_quality = 0; base_quality < 94; base_quality++) {
                for(map_quality = 0; map_quality < 256; map_quality++) {
                    double base_error = expPhred(base_quality);
                    double mapping_error = expPhred(map_quality);
                    double prob_if_mismapped = 0.25 * mapping_error;
                    probability_call_correct[base][called_base][base_quality][map_quality] = base == called_base ? (1.0 - base_error) : base_error / 3.0;
                    probability_call_correct[base][called_base][base_quality][map_quality] *= (1.0 - mapping_error);
                    probability_call_correct[base][called_base][base_quality][map_quality] += 0.25 * mapping_error;
                }
            }
        }
    }
}

#endif
