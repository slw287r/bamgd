# GC~Depth plot of microbial bam using cairographics

## Usage
```
./bamgd -i ../dat/573.bam -o ../dat/573.png -r ../dat/573.fa.gz
```
## Demo

Normal depth with dups             |  Normal depth dup-ed and dedup-ed via `-d` or `--dup`
:-------------------------:|:-------------------------:
<img src="https://github.com/user-attachments/assets/6fd8f188-2476-47b4-8305-8cab621c16cf" width=512></img>  |  <img src="https://github.com/user-attachments/assets/1ce990b1-8e53-4dcd-9256-b6685559fc1d" width=512></img>


Ultra-deepp depth             |  Ultra high depth logscaled via `-l` or `--log` 
:-------------------------:|:-------------------------:
<img src="https://github.com/user-attachments/assets/41e0a3aa-12a7-4779-914e-891fe32ff028" width=512></img> | <img src="https://github.com/user-attachments/assets/32437c90-4f37-40c5-9bd9-fee6c4d744a8" width=512></img>


| Ultra high depth logscaled with dup-ed and dedup-ed depth via `-ld` or `--log --dup` |
:-------------------------:
 |<img src="https://github.com/user-attachments/assets/3af7ef00-519a-4721-a322-edb2e348edd5" width=512></img>|
