/**
 * @author See Contributors.txt for code contributors and overview of BadgerDB.
 *
 * @section LICENSE
 * Copyright (c) 2012 Database Group, Computer Sciences Department, University of Wisconsin-Madison.
 */

#include "btree.h"
#include "filescan.h"
#include "exceptions/bad_index_info_exception.h"
#include "exceptions/bad_opcodes_exception.h"
#include "exceptions/bad_scanrange_exception.h"
#include "exceptions/no_such_key_found_exception.h"
#include "exceptions/scan_not_initialized_exception.h"
#include "exceptions/index_scan_completed_exception.h"
#include "exceptions/file_not_found_exception.h"
#include "exceptions/end_of_file_exception.h"


//#define DEBUG

namespace badgerdb
{

// -----------------------------------------------------------------------------
// BTreeIndex::BTreeIndex -- Constructor
// -----------------------------------------------------------------------------

BTreeIndex::BTreeIndex(const std::string & relationName,
		std::string & outIndexName,
		BufMgr *bufMgrIn,
		const int attrByteOffset,
		const Datatype attrType)
{	
	std::ostringstream idxStr;
	idxStr << relationName << '.' << attrByteOffset;
	std::string indexName = idxStr.str(); // indexName is the name of the index file.

	bool create_new = !file->exists(indexName);
	file = new BlobFile(indexName, create_new);
	outIndexName = indexName;
	bufMgr = bufMgrIn;
	scanExecuting = false;

	IndexMetaInfo *meta;
	LeafNodeInt *root;
	Page *headerPage, *rootPage;
	RecordId rid;

	if (create_new) {
		FileScan relation(relationName, bufMgr);

		bufMgr->allocPage(file, headerPageNum, headerPage);
		bufMgr->allocPage(file, rootPageNum, rootPage);

		meta = (IndexMetaInfo *)headerPage;
		strcpy(meta->relationName, relationName.c_str());
		meta->attrByteOffset = attrByteOffset;
		meta->attrType = attrType;
		meta->rootPageNo = rootPageNum;

		root = (LeafNodeInt *)rootPage;
		relation.scanNext(rid);
		root -> keyArray[0] = *((int*)relation.getRecord().c_str() + attrByteOffset);
		root -> ridArray[0] = rid;
		root->rightSibPageNo = 0;

		bufMgr->unPinPage(file, headerPageNum, true);
		bufMgr->unPinPage(file, rootPageNum, true);

		try {
			while (true) {
				relation.scanNext(rid);
				insertEntry(relation.getRecord().c_str() + attrByteOffset, rid);
			}
		} catch(EndOfFileException &e) {
			bufMgr->flushFile(file);
		}

	} else {
		bufMgr->readPage(file, 1, headerPage);
		
		meta = (IndexMetaInfo *)headerPage;
		rootPageNum = meta->rootPageNo;

		bufMgr->unPinPage(file, 1, false);
	}
}

// -----------------------------------------------------------------------------
// BTreeIndex::~BTreeIndex -- destructor
// -----------------------------------------------------------------------------

BTreeIndex::~BTreeIndex()
{	
	//update the file in the disk, 
	//delete the file in buffer manager
	bufMgr->flushFile(file);
	scanExecuting = false;
	delete file;
	file = nullptr;
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertEntry
// -----------------------------------------------------------------------------

void BTreeIndex::insertEntry(const void *key, const RecordId rid) 
{
	RIDKeyPair<int> leafPair;
	bool isLeaf = rootPageNum == 2 ? true : false;

	leafPair.set(rid, *((int *)key));
	insertHelper(leafPair, rootPageNum, isLeaf);
}

PageKeyPair<int> *BTreeIndex::insertHelper(RIDKeyPair<int> leafPair, PageId nodeId, bool isLeaf) {
	Page *newPage;
	PageKeyPair<int> *pushedPair;
	bufMgr->readPage(file, nodeId, newPage);

	if (isLeaf) {
		LeafNodeInt *Node = (LeafNodeInt *)newPage;
		if (Node->ridArray[INTARRAYLEAFSIZE-1].page_number == 0) { // have space
			insertLeaf(leafPair, Node);
			bufMgr -> unPinPage(file, nodeId, true);
			return nullptr;
		} else { // split
			pushedPair = splitLeaf(leafPair, nodeId);
			bufMgr -> unPinPage(file, nodeId, true);
			return pushedPair;
		}
	} else {
		NonLeafNodeInt *Node = (NonLeafNodeInt *)newPage;
		bool isChildLeaf = Node->level == 1 ? true : false;
		for (int i = 0; i < INTARRAYNONLEAFSIZE; i++) {
			if (i == 0 && leafPair.key < Node->keyArray[i]) {
				pushedPair = insertHelper(leafPair, Node->pageNoArray[0], isChildLeaf);
				break;
			}
			if (Node->keyArray[i-1] <= leafPair.key && leafPair.key < Node->keyArray[i]) {
				pushedPair = insertHelper(leafPair, Node->pageNoArray[i], isChildLeaf);
				break;
			} else if (Node->keyArray[i] == 0) {
				pushedPair = insertHelper(leafPair, Node->pageNoArray[i], isChildLeaf);
				break;
			} else if (i == INTARRAYNONLEAFSIZE-1 && Node->keyArray[i] <= leafPair.key) {
				pushedPair = insertHelper(leafPair, Node->pageNoArray[i+1], isChildLeaf);
				break;
			}
		}
		if (pushedPair != nullptr) {
			if (Node->pageNoArray[INTARRAYNONLEAFSIZE] == 0) {
				insertNonLeaf(*pushedPair, Node);
				bufMgr -> unPinPage(file, nodeId, true);
				return nullptr;
			} else {
				pushedPair = splitNonLeaf(pushedPair, nodeId);
				bufMgr -> unPinPage(file, nodeId, true);
				return pushedPair;
			}
		}
		bufMgr -> unPinPage(file, nodeId, true);
		return nullptr;
	}
}

void BTreeIndex::insertLeaf(RIDKeyPair<int> leafPair, LeafNodeInt *Node) {
	RIDKeyPair<int> cmpPair;
	for (int i = 0; i < INTARRAYLEAFSIZE; i++) {
		cmpPair.set(Node->ridArray[i], Node->keyArray[i]);

		if (leafPair < cmpPair) {
			Node->keyArray[i] = leafPair.key;
			Node->ridArray[i] = leafPair.rid;
			leafPair = cmpPair;
		} else if (Node->ridArray[i].page_number == 0) {
			Node->keyArray[i] = leafPair.key;
			Node->ridArray[i] = leafPair.rid;
			break;
		}
	}
}

void BTreeIndex::insertNonLeaf(PageKeyPair<int> nonLeafPair, NonLeafNodeInt *Node) {
	PageKeyPair<int> tmpPair;
	for (int i = 0; i < INTARRAYNONLEAFSIZE; i++) {
		if (nonLeafPair.key < Node->keyArray[i]) {
			tmpPair.set(Node->pageNoArray[i + 1], Node->keyArray[i]);
			Node->keyArray[i] = nonLeafPair.key;
			Node->pageNoArray[i + 1] = nonLeafPair.pageNo;
			nonLeafPair = tmpPair;
		} else if (Node->pageNoArray[i+1] == 0) {
			Node->keyArray[i] = nonLeafPair.key;
			Node->pageNoArray[i + 1] = nonLeafPair.pageNo;
			break;
		}
	}
}

PageKeyPair<int> *BTreeIndex::splitLeaf(RIDKeyPair<int> leafPair, PageId nodeId) {
	LeafNodeInt *Node;
	LeafNodeInt *siblingNode;
	PageId pageId;
	Page *newPage;

	bufMgr->readPage(file, nodeId, newPage);
	Node = (LeafNodeInt *)newPage;
	bufMgr->allocPage(file, pageId, newPage);
	siblingNode = (LeafNodeInt *)newPage;

	// bufMgr->unPinPage(file, nodeId, true);
	// bufMgr->unPinPage(file, pageId, true);

	int mid = INTARRAYLEAFSIZE / 2;

	if (Node->rightSibPageNo != 0) {
		siblingNode->rightSibPageNo = Node->rightSibPageNo;
	}
	Node->rightSibPageNo = pageId;

	for (int i = mid; i < INTARRAYLEAFSIZE; i++) {
		siblingNode->keyArray[i - mid] = Node->keyArray[i];
		siblingNode->ridArray[i - mid] = Node->ridArray[i];
		Node->keyArray[i] = 0;
		Node->ridArray[i].page_number = 0;
		Node->ridArray[i].slot_number = 0;
	}

	RIDKeyPair<int> cmpPair;
	cmpPair.set(siblingNode->ridArray[0], siblingNode->keyArray[0]);
	if (leafPair < cmpPair) {
		insertLeaf(leafPair, Node);
	} else {
		insertLeaf(leafPair, siblingNode);
	}

	bufMgr->unPinPage(file, nodeId, true);
	bufMgr->unPinPage(file, pageId, true);
	return pushUp(nodeId, Node->rightSibPageNo, siblingNode->keyArray[0]);
}

PageKeyPair<int> *BTreeIndex::splitNonLeaf(PageKeyPair<int> *nonLeafPair, PageId nodeId) {
	NonLeafNodeInt *Node;
	NonLeafNodeInt *siblingNode;
	PageId pageId;
	Page *newPage;
	int midKey;

	bufMgr->readPage(file, nodeId, newPage);
	Node = (NonLeafNodeInt *)newPage;
	bufMgr->allocPage(file, pageId, newPage);
	siblingNode = (NonLeafNodeInt *)newPage;

	// bufMgr->unPinPage(file, nodeId, true);
	// bufMgr->unPinPage(file, pageId, true);

	int mid = INTARRAYNONLEAFSIZE / 2;

	siblingNode->level = Node->level;
	// if (Node->keyArray[mid-1] < nonLeafPair->key && nonLeafPair->key < Node->keyArray[mid]) {

	// }
	siblingNode->pageNoArray[0] = Node->pageNoArray[mid];
	for (int i = mid; i < INTARRAYNONLEAFSIZE; i++) {
		siblingNode->keyArray[i - mid] = Node->keyArray[i];
		siblingNode->pageNoArray[i - mid + 1] = Node->pageNoArray[i + 1];
		Node->keyArray[i] = 0;
		Node->pageNoArray[i + 1] = 0;
	}
	// Node->pageNoArray[mid] = nonLeafPair->pageNo;

	if (siblingNode->keyArray[0] > nonLeafPair->key) {
		insertNonLeaf(*nonLeafPair, Node);
		midKey = Node->keyArray[mid];
		Node->keyArray[mid] = 0;
		Node->pageNoArray[mid+1] = 0;
	} else {
		insertNonLeaf(*nonLeafPair, siblingNode);
		midKey = Node->keyArray[mid-1];
		Node->keyArray[mid-1] = 0;
		Node->pageNoArray[mid] = 0;
	}

	bufMgr->unPinPage(file, nodeId, true);
	bufMgr->unPinPage(file, pageId, true);
	// PageKeyPair<int> *rightPair = new PageKeyPair<int>; // push up
	// rightPair->set(pageId, midKey);
	return pushUp(nodeId, pageId, midKey); // push up
}

PageKeyPair<int> *BTreeIndex::pushUp(PageId leftNodeId, PageId rightNodeId, int key) {
	NonLeafNodeInt *root;
	Page *newPage;

	if (leftNodeId == rootPageNum) {  // new root node
		PageId rootId;
		IndexMetaInfo *meta;
		bufMgr->allocPage(file, rootId, newPage);
		root = (NonLeafNodeInt *)newPage;
		root->level = rootPageNum == 2 ? 1 : 0;
		root->keyArray[0] = key;
		root->pageNoArray[0] = leftNodeId;
		root->pageNoArray[1] = rightNodeId;

		rootPageNum = rootId; // modify rootPageNum
		bufMgr->readPage(file, headerPageNum, newPage);
		meta = (IndexMetaInfo *)newPage;
		meta->rootPageNo = rootId;
		bufMgr->unPinPage(file, rootId, true);
		bufMgr->unPinPage(file, headerPageNum, true);
		return nullptr;
	} else {
		PageKeyPair<int> *rightPair = new PageKeyPair<int>;
		rightPair->set(rightNodeId, key);
		return rightPair;
	}
}

// -----------------------------------------------------------------------------
// BTreeIndex::startScan
// -----------------------------------------------------------------------------

void BTreeIndex::startScan(const void* lowValParm,
				   const Operator lowOpParm,
				   const void* highValParm,
				   const Operator highOpParm)
{
	this->scanExecuting =false;
	this->highValInt = *((int *)highValParm);
	this->lowValInt = *((int *)lowValParm);
	this->lowOp = lowOpParm;
	this->highOp = highOpParm;
	if ((lowOpParm != GT && lowOpParm != GTE ) 
	|| (highOpParm != LT && highOpParm != LTE))
	{
		throw BadOpcodesException();
	}else if (this->lowValInt>this->highValInt)
	{
		throw BadScanrangeException(); 
	}

	NonLeafNodeInt *curNode = NULL; // current node 
	this->currentPageNum= rootPageNum; // current page ID 
	int arri;
	int layer = 0;
	LeafNodeInt *leafNode = NULL;
	if (rootPageNum == 2)
	{
		layer = 1;
		bufMgr->readPage(file,this->currentPageNum,(Page *&)leafNode); //will pin the page
		bufMgr->unPinPage(file, this->currentPageNum,false); 
	} else {
		bufMgr->readPage(file, this->currentPageNum,(Page *&)curNode);
		bufMgr->unPinPage(file,this->currentPageNum,false); //unpin the page
	}
	while(layer==0) // when level equal 1, we reach level that is above leaf layer
	{
		for (arri=0; arri<INTARRAYNONLEAFSIZE;arri++){
			if (curNode->keyArray[arri]== 0 ||	// 0 is when there is nothing in that element.
			curNode->keyArray[arri] > this->lowValInt){
				break;
			}
		}
		this->currentPageNum = curNode->pageNoArray[arri];
		layer= curNode->level;
		if (layer==0){
			bufMgr->readPage(file,this->currentPageNum,(Page *&)curNode); //will pin the page
			bufMgr->unPinPage(file, this->currentPageNum,false); 
		}else{
			bufMgr->readPage(file,this->currentPageNum,(Page *&)leafNode); //will pin the page
			bufMgr->unPinPage(file, this->currentPageNum,false); 
		}
		
	}

	//in leaf node
	for (arri=0; arri<INTARRAYLEAFSIZE-1;arri++){
		if (leafNode->keyArray[arri] == 0 && leafNode->keyArray[arri+1] == 0){	// 0 is when there is nothing in that element.
			throw NoSuchKeyFoundException();
		} else {	// 0 is when there is nothing in that element.
			if (this->lowOp == GT){
				if (leafNode->keyArray[arri] > this->lowValInt)
					break;
			}else if (this->lowOp == GTE){
				if (leafNode->keyArray[arri] >= this->lowValInt)
					break;
			}
		}
	}
	if (arri == INTARRAYLEAFSIZE-1 && leafNode->keyArray[arri] == 0) {
		throw NoSuchKeyFoundException();
	}
	//In this step we get a key must higher than lower bound
	//now checking this key is lower than higher bound
	switch (this->highOp)
	{
	case LT:
		if (leafNode->keyArray[arri] > this->highValInt)
			throw NoSuchKeyFoundException();
		break;
	case LTE:
		if (leafNode->keyArray[arri] >= this->highValInt)
			throw NoSuchKeyFoundException();
		break;
	default:
		break;
	}
	this->scanExecuting = true;
	this->currentPageData = (Page *)leafNode;
	this->nextEntry = arri;
}

// -----------------------------------------------------------------------------
// BTreeIndex::scanNext
// -----------------------------------------------------------------------------

void BTreeIndex::scanNext(RecordId& outRid) 
{
	if (this->scanExecuting==false)
		throw ScanNotInitializedException();
	LeafNodeInt *node = (LeafNodeInt * )this->currentPageData;
	//if current entry is invalid 
	if (this->nextEntry>= INTARRAYLEAFSIZE || 
		(this->nextEntry < INTARRAYLEAFSIZE-1 && node->keyArray[this->nextEntry] == 0 && node->keyArray[this->nextEntry+1] == 0) || 
		(this->nextEntry == INTARRAYLEAFSIZE-1 && node->keyArray[this->nextEntry] == 0)){	// true need to switch page, 0 is when there is nothing in that element.
		if(node->rightSibPageNo == 0) //or try this (PageId)EMPTY_SLOT
			throw IndexScanCompletedException();
		//now switch to next page
		bufMgr->readPage(file, node->rightSibPageNo, (Page *&)this->currentPageData);
		bufMgr->unPinPage(file, node->rightSibPageNo,false); 
		this->nextEntry = 0;
		this->currentPageNum = node->rightSibPageNo;
		node = (LeafNodeInt * )this->currentPageData;
	}
	//In this step we get a key must higher than lower bound
	//now checking this key is lower than higher bound
	switch (this->highOp)
	{
	case LT:
		if (node->keyArray[this->nextEntry] >= this->highValInt)
			throw IndexScanCompletedException();
		break;
	case LTE:
		if (node->keyArray[this->nextEntry] > this->highValInt)
			throw IndexScanCompletedException();
		break;
	default:
		break;
	}
	outRid = node->ridArray[this->nextEntry];
	//move to next entry
	this->nextEntry +=1;
}

// -----------------------------------------------------------------------------
// BTreeIndex::endScan
// -----------------------------------------------------------------------------
//
void BTreeIndex::endScan() 
{
	//check is scannig status
	if(this->scanExecuting==false)
		throw ScanNotInitializedException();

	//clear the object
	this->nextEntry = 0;
	this->currentPageData = nullptr;
	this->currentPageNum = 0;
	this->scanExecuting = false;
}

}