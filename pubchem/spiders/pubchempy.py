import scrapy
import json
# from pubchem.items import PubchemItem
from scrapy_splash import SplashRequest
import time
import sys
import re


class PubchempySpider(scrapy.Spider):
    http_user = 'ad4d3c8537874b8fb4104412c68845ba' 
    name = 'pubchempy'
    allowed_domains = ['pubchem.ncbi.nlm.nih.gov']
    cid_list = []
    # descrição do script na linguagem LUA para que seja aberta a página para poder extrair dados inicializados em javascript
    script = '''
        function main(splash, args)
            splash.private_mode_enabled = false
            url = args.url
            assert(splash:go(args.url))
            assert(splash:wait(8))
            splash:set_viewport_full()
            return {
                html = splash:html()
            }
        end
    '''
    cas_list_obj = cas_list = ['75-07-0']

    def start_requests(self, cas_list_obj=cas_list_obj):
        host = 'https://pubchem.ncbi.nlm.nih.gov/#query='
        start_urls = []
        # início das requisições para serem feitas utilizando uma lista de CAS ou inchikeys
        # essas requisições abrem páginas de pesquisa utilizando os dados da lista no site do pubchem
        cas_list = cas_list_obj

        for each_cas in cas_list:
            url_query = host + each_cas
            start_urls.append(url_query)
        for each_host in start_urls:
            yield SplashRequest(url=each_host, callback=self.get_cid_or_sid, endpoint="execute", args={
                'lua_source': self.script
            })

    def get_cid_or_sid(self, response):
        contador = 0
        # item = PubchemItem()
        # verifica se a linha com o CID contém também SID
        compound_cid = response.xpath(
            "//div[contains(@class, 'flex-grow-1 p-md-left')]/div[2]//div[2]//a//span/text()")
        cas_desejado = response.xpath(
            "//span[contains(@class, 'highlight')]/text()").get()

        if len(compound_cid) > 0:
            for compound in response.xpath("//div[contains(@class, 'flex-grow-1 p-md-left')]/div[2]"):
                contador += 1
                if contador >= 2:
                    pass
                else:
                    # compound_cid = compound.xpath(".//div[2]//a//span/text()").get()
                    compound_cid = compound.xpath(".//a//span/text()").get()
        else:
            for compound in response.xpath("//div[contains(@class, 'flex-grow-1 p-md-left')]/div[2]"):
                contador += 1
                if contador >= 2:
                    pass
                else:
                    # compound_cid = compound.xpath(".//div[2]//a//span/text()").get()
                    compound_cid = compound.xpath(".//a//span/text()").get()

        yield {
            'cas_desejado': cas_desejado,
            'compound_cid': compound_cid
        }

'''
        # checa se é SUBSTANCE SID ou COMPOUND CID
        name_id = response.xpath(
            "//div[contains(@class, 'flex-grow-1 p-md-left')]/div[2]")
        count = 0
        for each_name_id in name_id:
            count += 1
            if count >= 2:
                pass
            else:
                name_id = each_name_id.xpath(".//span/text()").get()
        # início do loop para extrair as informações do pubchem a partir do CAS
        if "CID" in name_id:
            host = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/'
            last_host = '/JSON/'
            json_url = host + str(compound_cid) + last_host

            yield response.follow(
                url=json_url,
                callback=self.get_pubchem_data,
                meta={'cas_desejado': cas_desejado},
                dont_filter=True)

        elif "SID" in name_id:
            host = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/substance/'
            last_host = '/JSON/'
            json_url = host + str(compound_cid) + last_host

            yield response.follow(
                url=json_url,
                callback=self.get_pubchem_data,
                meta={'cas_desejado': cas_desejado},
                dont_filter=True)

    def unique(self, sequence):
        seen = set()
        return [x for x in sequence if not (x in seen or seen.add(x))]

# ---------------------------------- métodos para a extração de dados a partir do json das substâncias/compostos
    def get_title(self, json_file):
        try:
            if json_file['Record']['RecordType'] == 'CID':
                title = json_file['Record']['RecordTitle']
            elif json_file['Record']['RecordType'] == 'SID':
                title = json_file['Record']['RecordTitle']
            return title
        except:
            title = 'erro na obtenção da title'
            return title

    def get_cas(self, json_file, cas_desejado):
        cas_pattern = re.compile(r'\b[1-9]{1}[0-9]{1,10}-\d{2}-\d\b|\bCAS-[1-9]{1}[0-9]{1,10}-\d{2}-\d\b')
        Lista_vazia_de_cas = []
        lista_vazia = []
        lista_de_cas_possiveis = []
        try:
            for each_node in json_file['Record']['Section']:
                if each_node['TOCHeading'] == "Names and Identifiers":
                    for r in each_node['Section']:
                        if r['TOCHeading'] == "Synonyms":
                            for s in r['Section']:
                                if s['TOCHeading'] == "Depositor-Supplied Synonyms":
                                    for w in s['Information']:
                                        for t in w['Value'].values():
                                            Lista_vazia_de_cas.append(t)
                                        for i in Lista_vazia_de_cas:
                                            for r in i:
                                                lista_vazia.append(r['String'])
            lista_de_cas_sem_repeticao = self.unique(lista_vazia)
            for i in lista_de_cas_sem_repeticao:
                if bool(re.match(cas_pattern, i)) is True:
                    # cas = i
                    # break
                    lista_de_cas_possiveis.append(i)
            if cas_desejado in lista_de_cas_possiveis:
                cas = cas_desejado
            return cas
        except:
            try:
                for each_node in json_file['Record']['Section']:
                    if each_node['TOCHeading'] == "Names and Identifiers":
                        for r in each_node['Section']:
                            if r['TOCHeading'] == "Other Identifiers":
                                for s in r['Section']:
                                    if s['TOCHeading'] == "CAS":
                                        for w in s['Information']:
                                            for t in w['Value'].values():
                                                Lista_vazia_de_cas.append(t)
                                            for i in Lista_vazia_de_cas:
                                                for r in i:
                                                    lista_vazia.append(
                                                        r['String'])
                lista_de_cas_sem_repeticao = self.unique(lista_vazia)
                # for i in lista_de_cas_sem_repeticao:
                #     if bool(re.match(cas_pattern, i)) is True:
                #         cas = i
                #         break
                # return cas
                for i in lista_de_cas_sem_repeticao:
                    if bool(re.match(cas_pattern, i)) is True:
                        # cas = i
                        # break
                        lista_de_cas_possiveis.append(i)
                if cas_desejado in lista_de_cas_possiveis:
                    cas = cas_desejado
                return cas
            except:
                cas = 'end'
                return cas

    def get_cas_from_record0(self, json_file, cas_desejado):
        cas_pattern = re.compile(r'\b[1-9]{1}[0-9]{1,10}-\d{2}-\d\b|\bCAS-[1-9]{1}[0-9]{1,10}-\d{2}-\d\b')
        lista_cas1 = []
        lista_vazia_cas1 = []
        lista_de_cas_possiveis = []
        try:
            for each_node in json_file['Record']['Section']:
                if each_node['TOCHeading'] == "Names and Identifiers":
                    for r in each_node['Section']:
                        if r['TOCHeading'] == "Other Identifiers":
                            for s in r['Section']:
                                if s['TOCHeading'] == "CAS":
                                    for w in s['Information']:
                                        for t in w['Value'].values():
                                            lista_cas1.append(t)
                                        for i in lista_cas1:
                                            for r in i:
                                                lista_vazia_cas1.append(
                                                    r['String'])
            lista_de_cas_sem_repeticao = self.unique(lista_vazia_cas1)
            # for i in lista_de_cas_sem_repeticao:
            #     if bool(re.match(cas_pattern, i)) is True:
            #         cas = i
            #         break
            # return cas
            for i in lista_de_cas_sem_repeticao:
                if bool(re.match(cas_pattern, i)) is True:
                    # cas = i
                    # break
                    lista_de_cas_possiveis.append(i)
            if cas_desejado in lista_de_cas_possiveis:
                cas = cas_desejado
            return cas
        except:
            cas = 'end'
            return cas

    # função para capturar CAS dos óleos minerais e naftas
    def get_cas_from_record1(self, json_file, cas_desejado):
        cas_pattern = re.compile(r'\b[1-9]{1}[0-9]{1,10}-\d{2}-\d\b|\bCAS-[1-9]{1}[0-9]{1,10}-\d{2}-\d\b')
        lista_cas1 = []
        lista_vazia_cas1 = []
        lista_de_cas_possiveis = []
        try:
            for each_node in json_file['Record']['Section']:
                if each_node['TOCHeading'] == "Identity":
                    for r in each_node['Section']:
                        if r['TOCHeading'] == "Depositor-Supplied Synonyms":
                            for w in r['Information']:
                                for t in w['Value'].values():
                                    lista_cas1.append(t)
                                for i in lista_cas1:
                                    for r in i:
                                        lista_vazia_cas1.append(r['String'])
            lista_de_cas_sem_repeticao = self.unique(lista_vazia_cas1)
            # for i in lista_de_cas_sem_repeticao[0:]:
            #     if bool(re.match(cas_pattern, i)) is True:
            #         cas = i
            #         break
            # return cas
            for i in lista_de_cas_sem_repeticao:
                if bool(re.match(cas_pattern, i)) is True:
                    # cas = i
                    # break
                    lista_de_cas_possiveis.append(i)
            if cas_desejado in lista_de_cas_possiveis:
                cas = cas_desejado
            return cas
        except:
            cas = 'end'
            return cas

    # função para capturar o nome dos óleos minerais e naftas
    def get_name_from_record1(self, json_file):
        # cas_pattern = re.compile(r'\b[1-9]{1}[0-9]{1,10}-\d{2}-\d\b')
        cas_pattern = re.compile(r'\b[1-9]{1}[0-9]{1,10}-\d{2}-\d\b|\bCAS-[1-9]{1}[0-9]{1,10}-\d{2}-\d\b')
        lista_nomes = []
        lista_nomes_vazia = []
        lista_de_cas_possiveis = []
        try:
            for each_node in json_file['Record']['Section']:
                if each_node['TOCHeading'] == "Identity":
                    for r in each_node['Section']:
                        if r['TOCHeading'] == "Depositor-Supplied Synonyms":
                            for w in r['Information']:
                                for t in w['Value'].values():
                                    lista_nomes.append(t)
                                for i in lista_nomes:
                                    for r in i:
                                        lista_nomes_vazia.append(r['String'])
            lista_de_nomes_sem_repeticao = self.unique(lista_nomes_vazia)
            for i in lista_de_nomes_sem_repeticao[0:]:
                if bool(re.match(cas_pattern, i)) is False:
                    nome = i
                    break
            return nome
            #--------------------------
            # for i in lista_de_cas_sem_repeticao:
            #     if bool(re.match(cas_pattern, i)) is True:
            #         # cas = i
            #         # break
            #         lista_de_cas_possiveis.append(i)
            # if cas_desejado in lista_de_cas_possiveis:
            #     cas = cas_desejado
            # return cas
        except:
            nome = 'end'
            return nome

    # função para capturar o tipo: se é óleo minerai ou nafta
    def get_type_from_record1(self, json_file):
        lista_nomes = []
        lista_nomes_vazia = []
        try:
            for each_node in json_file['Record']['Section']:
                if each_node['TOCHeading'] == "Identity":
                    for r in each_node['Section']:
                        if r['TOCHeading'] == "Depositor-Supplied Synonyms":
                            for w in r['Information']:
                                for t in w['Value'].values():
                                    lista_nomes.append(t)
                                for i in lista_nomes:
                                    for r in i:
                                        lista_nomes_vazia.append(r['String'])
            lista_de_nomes_sem_repeticao = self.unique(lista_nomes_vazia)
            for i in lista_de_nomes_sem_repeticao:
                if 'diesel' in i.lower():
                    tipo = 'diesel'
                elif 'oil' and 'mineral' in i.lower():
                    tipo = 'oleo_mineral'
                elif 'coal' and 'mineral' in i.lower():
                    tipo = 'carvao_mineral'
                elif 'bituminous' in i.lower():
                    tipo = 'carvao_mineral'
                elif 'lignite' in i.lower():
                    tipo = 'carvao_mineral'
                elif 'naphtha' in i.lower():
                    tipo = 'derivados_de_petroleo'
                elif 'naphthenic' in i.lower():
                    tipo = 'derivados_de_petroleo'
                elif 'bitume' in i.lower():
                    tipo = 'derivados_de_petroleo'
                elif 'fuel' in i.lower():
                    tipo = 'derivados_de_petroleo'
                elif 'kerosene' in i.lower():
                    tipo = 'derivados_de_petroleo'
                elif 'kerosine' in i.lower():
                    tipo = 'derivados_de_petroleo'
                elif 'petroleum' in i.lower():
                    tipo = 'derivados_de_petroleo'
                elif 'distillates' in i.lower():
                    tipo = 'derivados_de_petroleo'
                elif 'solvent' in i.lower():
                    tipo = 'derivados_de_petroleo'
                elif 'paraffin' in i.lower():
                    tipo = 'derivados_de_petroleo'
            return tipo
        except:
            tipo = 'end'
            return tipo

    def get_tipo(self, json_file):
        try:
            if json_file['Record']['RecordType'] == 'CID':
                cid = 'compound'
            elif json_file['Record']['RecordType'] == 'SID':
                cid = 'substance'
            return cid
        except:
            cid = 'erro na obtenção da cid'
            return cid

    def get_cid_on_json(self, json_file):
        try:
            if json_file['Record']['RecordType'] == 'CID':
                cid = json_file['Record']['RecordNumber']
            elif json_file['Record']['RecordType'] == 'SID':
                cid = json_file['Record']['RecordNumber']
            return cid
        except:
            cid = 'erro na obtenção da cid'
            return cid

    def get_iupac_name(self, json_file):
        void = []
        iupac_list = []
        try:
            for each_node in json_file['Record']['Section']:
                if each_node['TOCHeading'] == "Names and Identifiers":
                    for r in each_node['Section']:
                        if r['TOCHeading'] == "Computed Descriptors":
                            for s in r['Section']:
                                if s['TOCHeading'] == "IUPAC Name":
                                    for w in s['Information']:
                                        for t in w['Value'].values():
                                            void.append(t)
                                        for i in void:
                                            for r in i:
                                                iupac_list.append(r['String'])
            iupac_name = str(iupac_list[0])
            return iupac_name
        except:
            iupac_name = 'error getting iupac name'
            return iupac_name

    def get_sinon(self, json_file):
        trash_pattern = re.compile(
            r'(\bApp\b|\bUnii\b|\bAcmc\b|\bFluka\b|\b[1-9]{1}[0-9]{1,10}-\d{2}-\d\b|[0-9]{2,}?[a-zA-Z]{1,})|([a-zA-Z_]{1,}?([0-9]{1,}))|(\#\d+|\#\w+|\w+\#|\d+\#)|([a-zA-Z]{1,}( |:|-|_)?[a-zA-Z]{1,}( |:|-|_)?([0-9]{2,}))', re.IGNORECASE)
        sinon_list = []
        sinon_list_strings = []
        melhor_lista = []
        contador = 0
        try:
            for each_node in json_file['Record']['Section']:
                if each_node['TOCHeading'] == "Names and Identifiers":
                    for r in each_node['Section']:
                        if r['TOCHeading'] == "Synonyms":
                            for s in r['Section']:
                                if s['TOCHeading'] == "Depositor-Supplied Synonyms":
                                    for w in s['Information']:
                                        for t in w['Value'].values():
                                            sinon_list.append(t)
                                        for i in sinon_list:
                                            for r in i:
                                                sinon_list_strings.append(
                                                    r['String'])
            for cada_sinonimo in sinon_list_strings:
                if contador < 5:
                    melhor_lista.append(cada_sinonimo)
                else:
                    pass
                contador += 1
            lista_de_cas_sem_repeticao = self.unique(melhor_lista)
            filtered = [
                i for i in lista_de_cas_sem_repeticao if not trash_pattern.match(i)]
            sinonimos = filtered
            return sinonimos
        except:
            sinonimos = 'error getting sinom'
            return sinonimos

    def get_inchikey(self, json_file):
        void = []
        inchikey_list = []
        try:
            for each_node in json_file['Record']['Section']:
                if each_node['TOCHeading'] == "Names and Identifiers":
                    for r in each_node['Section']:
                        if r['TOCHeading'] == "Computed Descriptors":
                            for s in r['Section']:
                                if s['TOCHeading'] == "InChI Key":
                                    for w in s['Information']:
                                        for t in w['Value'].values():
                                            void.append(t)
                                        for i in void:
                                            for r in i:
                                                inchikey_list.append(
                                                    r['String'])
            inchikey = str(inchikey_list[0])
            return inchikey
        except:
            inchikey = 'error getting inchikey name'
            return inchikey

# --------------------------------- fim dos métodos usados para extração de dados a partir do json das substâncias/compostos

    def get_pubchem_data(self, response, cas_list_obj=cas_list_obj):
        cas_list = cas_list_obj
        cas_desejado = response.meta.get('cas_desejado')

        host = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/'
        json_resp = json.loads(response.body)
        last_host = str(self.get_tipo(json_resp)) + str('/') + \
            str(self.get_cid_on_json(json_resp)) + str('/json')
        cas1 = self.get_cas(json_resp, cas_desejado)
        cas2 = self.get_cas_from_record0(json_resp, cas_desejado)
        cas3 = self.get_cas_from_record1(json_resp, cas_desejado)
        if cas1 in cas_list:
            print("executando a parte do cas1")
            yield {
                'Tipo': self.get_tipo(json_resp),
                'URL do agente quimico': host + last_host,
                'CAS': self.get_cas(json_resp, cas_desejado),
                'Nome no cabecalho': self.get_title(json_resp),
                'Nome IUPAC': self.get_iupac_name(json_resp),
                'Sinonimos': self.get_sinon(json_resp),
                'inchikey': self.get_inchikey(json_resp),
                'Full json': json_resp
            }
        elif cas2 in cas_list:
            print("executando a parte do cas2")
            yield {
                'Tipo': self.get_tipo(json_resp),
                'URL do agente quimico': host + last_host,
                'CAS': self.get_cas_from_record0(json_resp, cas_desejado),
                'Nome no cabecalho': self.get_title(json_resp),
                'Nome IUPAC': self.get_iupac_name(json_resp),
                'Sinonimos': self.get_sinon(json_resp),
                'inchikey': self.get_inchikey(json_resp),
                'Full json': json_resp
            }
        elif cas3 in cas_list:
            print("executando a parte do cas3")
            yield {
                'Tipo': self.get_tipo(json_resp),
                'URL do agente quimico': host + last_host,
                'CAS': self.get_cas_from_record1(json_resp, cas_desejado),
                'Nome no cabecalho': self.get_name_from_record1(json_resp),
                'Tipo do refino': self.get_type_from_record1(json_resp),
                'Full json': json_resp
            }
        else:
            pass
'''