var xhr = new XMLHttpRequest();
var db;
var syntenic_data = {}

function testWorker(){
    var worker = new Worker("worker.js");

    worker.postMessage({
        url: 'insr.sqlite',
        sql: "SELECT DISTINCT gene_family_id FROM gene_family"
    });

    worker.onmessage = function(event) {
        worker.terminate();
    };

    worker.onerror = function(e) {console.log("Worker error: ", e)};
}


function setDB(db_name, callback){
    xhr.open('GET', db_name, true);
    xhr.responseType = 'arraybuffer';
    xhr.onload = function(e) {
        var uInt8Array = new Uint8Array(this.response);
        db = new SQL.Database(uInt8Array);
        // var gene_family_id = db.exec("SELECT DISTINCT gene_family_id FROM gene_family");

        // var html = "<SELECT id='gene_family' onchange=get_Genes_for_family()>"
        // for(var i=0; i<gene_family_id[0].values.length; i++){
        //     html += "<option value="+gene_family_id[0].values[i]+">"+gene_family_id[0].values[i]+"</option>"
        // }
        // html += "</SELECT>"

        // document.getElementById("families").innerHTML = html

        return callback()
    };
    xhr.send();
}



function get_Genes_for_family(gene_family_id, protein_id, gene_id){
    syntenic_data = {};

    if(!gene_family_id)
    {
     gene_family_id = 1;   
    }

    var gene_tree = db.exec("SELECT gene_family_member.protein_id, protein_alignment, transcript.gene_id, gene_json, transcript.protein_sequence FROM gene_family_member JOIN transcript ON gene_family_member.protein_id=transcript.protein_id JOIN gene ON transcript.gene_id=gene.gene_id WHERE gene_family_id= "+gene_family_id);
    var genes= {}
    var cigars= {}
    var sequence= {}

    for(var i=0; i<gene_tree[0].values.length; i++){
        cigars[gene_tree[0].values[i][0]] = gene_tree[0].values[i][1]
        genes[gene_tree[0].values[i][2]] = JSON.parse(gene_tree[0].values[i][3])
        sequence[gene_tree[0].values[i][0]] = gene_tree[0].values[i][4]
    }

    syntenic_data["member"] = genes
    syntenic_data["cigar"] = cigars
    syntenic_data["sequence"] = sequence

    syntenic_data["ref"] = gene_id ? gene_id : gene_tree[0].values[0][2]
    
    syntenic_data["protein_id"] = protein_id ? protein_id : gene_tree[0].values[0][0]

    syntenic_data["tree"] = get_GeneTree_for_family(gene_family_id)

    document.getElementById("gene_tree").innerHTML = "";
    removePopup()
    init(syntenic_data, "#settings_div", "#filter" , "#sliderfilter")


    drawTree(syntenic_data.tree, "#gene_tree", newpopup)
    
}

function get_Genes_by_ID(gene_id){
    var gene = db.exec("SELECT gene_json FROM gene_information where gene_id = '"+gene_id+"'");

    var gene_json
    for(var i=0; i<gene[0].values.length; i++){
        gene_json = JSON.parse(gene[0].values[i]);
    }
    return gene_json;
}



function get_CIGAR_for_gene(gene_id){

    var e = document.getElementById("gene_family");
    var gene_family_id = e.options[e.selectedIndex].value;

    var gene_tree = db.exec("SELECT protein_id, protein_alignment FROM gene_family where gene_family_id = "+gene_family_id+" and gene_id='"+gene_id+"'");

    var gene_cigar = {}
   
    for(var i=0; i<gene_tree[0].values.length; i++){
        gene_cigar[gene_tree[0].values[i][0]] = gene_tree[0].values[i][1]
    }

    return gene_cigar;
}

function search_for_gene(){

    var e = document.getElementById("seach_box");
    var keyword = e.value;

    var gene_id = db.exec("SELECT t.protein_id, g.species, g.gene_symbol, t.gene_id FROM gene g, transcript t where g.gene_symbol LIKE '%"+keyword+"%' and g.gene_id = t.gene_id and t.protein_id not null");

    var results = ""

    console.log(gene_id)
    for(var i=0; i<gene_id[0].values.length; i++){
        results += "<div style='margin: 5px; background: lightgray none repeat scroll 0% 0%; cursor: pointer; padding: 5px;' onclick=get_GeneFamily_by_Gene('"+gene_id[0].values[i][0]+"','"+gene_id[0].values[i][3]+"')> Protein:"+gene_id[0].values[i][0]+" <br>Species: " +gene_id[0].values[i][1]+" <br>Gene: " +gene_id[0].values[i][2]+" (" +gene_id[0].values[i][3]+")</div>"
    }

    document.getElementById("seach_results").innerHTML = results;
}

function get_GeneTree_for_family(gene_family_id){

    var gene_tree = db.exec("SELECT gene_tree FROM gene_family where gene_family_id = "+gene_family_id);

    var html = ""
    for(var i=0; i<gene_tree[0].values.length; i++){
        html += gene_tree[0].values[i]
    }

    return html;

}

function get_GeneFamily_by_Gene(protein_id, gene_id){
    var gene_family = db.exec("SELECT gene_family_id FROM gene_family_member where protein_id = '"+protein_id+"'");
    get_Genes_for_family(gene_family[0].values[0][0], protein_id, gene_id)
}
