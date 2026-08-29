/* crc32.h : CRC-32 (IEEE 802.3, the zlib/PNG polynomial) for the .gkmk payload checksum. */
#pragma once
#include <cstddef>
#include <cstdint>

inline uint32_t gkm_crc32(const void *data, size_t n, uint32_t crc = 0) {
	static uint32_t table[256];
	static bool init = false;
	if (!init) {
		for (uint32_t i = 0; i < 256; i++) {
			uint32_t c = i;
			for (int k = 0; k < 8; k++) c = (c & 1) ? (0xEDB88320u ^ (c >> 1)) : (c >> 1);
			table[i] = c;
		}
		init = true;
	}
	const unsigned char *p = (const unsigned char *)data;
	crc = ~crc;
	for (size_t i = 0; i < n; i++) crc = table[(crc ^ p[i]) & 0xFF] ^ (crc >> 8);
	return ~crc;
}
