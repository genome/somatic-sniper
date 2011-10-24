#pragma once

#include "output_format.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void output_classic_header(FILE *fh, const header_data_t *h);
void output_classic(FILE *fh, const sniper_output_t *p);

#ifdef __cplusplus
}
#endif /* __cplusplus */
