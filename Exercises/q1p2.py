C = "glnonppfnasjdfzfoixaupzudoigrlnjaoslrjufglrdgatpdggnjrtglnzafpzaobuazagjnzgglntprrblabvabnglnufobuazonajpxsapvaobenfoierdobbruoglnjfynjglnropxglfoitrjfguazgrsrvngraobuafgtrjglngdjortglngfbnglnznajnaslrtglnglavnzzgjngslnbentrjndzpfqnglnenifoofoirtaofognjvfoaepnuagnjuaxfoglnrttfoiglnznaaobglnzqxunjnunpbnbgringlnjufglrdgahrfogaobfoglnpdvfordzzkasnglngaoonbzafpzrtglneajinzbjftgfoidkufglglngfbnznnvnbgrzgaobzgfppfojnbspdzgnjzrtsaoyazzlajkpxknaqnbufglipnavzrtyajofzlnbzkjfgzalamnjnzgnbroglnpruzlrjnzglagjaordggrznafoyaofzlfoitpagonzzglnafjuazbajqaerynijaynznobaobtajglnjeasqzgfppznnvnbsrobnoznbfogravrdjotdpiprrvejrrbfoivrgfropnzzrynjglnefiinzgaobglnijnagnzggruoronajgl"

C = "xmeagiyiyzaglyagiknnyxciazeykrykrkfyktpxaoiacikorgpmazxiiacyahkragjixolhaznriykikfyxnrfkoiyxoqkorpaccxjneraxipzkfixfknneeaghagnroifaocikoinezgokfzaccmanqciarkehyafnkxtiykikfyxnrraoiqoahkoeiyxolkfyxnrcjzkxocikzicmgofixaoxolkijxziykorykcktaolcixictkoexomkoifaowangixaociyagckorcamraztkoikiatcxoiahyxfylarykcpgiktecixfpaccxjxnxiemazoaixfxolkokrgnickfikormxlgzxolagixicpgzpazigpiakjagixicpzxtkzecfyaanrkeckfyxnriyxoqcokigzknneaoneampnkejgitkoekmaztampnkefaoikxocrxcfxpnxokzemkfiazceagfkoiraiyxcaziykipgiceagagicyahckfyxnriykixitgciiyxoqpzkfixfknneazmkxnoahxmiyzaglyagifyxnryaarkjzkxoykcoaappacxixaoxixcpnkxoiykixihxnnkiikxokpacxixaoamcikigcugakchxiyagzazrxokzekoxtknctkoqoahcoaihyekfahralaznxaohkcoaijazohxiykjzkxoaokpkzhxiyagzchyecgfykoxtkncfkooaikrrcgjizkfiazajikxomzatjaaqckorcfyaanxoliykipkzktagoipacxixaohyxfytkoyanrciarke"

character_counts = {}
for c in C:
    try:
        character_counts[c] += 1
    except KeyError:
        character_counts[c] = 1
character_counts = {k: v for k, v in sorted(character_counts.items(), key=lambda item: item[1], reverse=True)}
print(character_counts)
print()

# Create a codebook by assigning the most common character to e and so on
letters = "eariotnslcudpmhgbfywkvxzjq"
# letters = "etaoinshrdlcumwfgypbvkjxqz"
# letters = "etainoshrdlucmfwygpbvkqjxz"

codebook = {}
i = 0
for key in character_counts:
    codebook[key] = letters[i]
    i += 1

print(codebook)
print()



M = ""
for c in C:
    M += codebook[c]

print(M)
print()



